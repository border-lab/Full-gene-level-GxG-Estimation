# Function_input.py  --  ANNOTATED READING COPY
#
# This is the same code as Function_input.py, line for line; only comments were
# added.  Comments marked "NOTE:" are review findings (stale text, edge cases,
# things worth knowing); none of them changes what the code computes.
#
# ============================================================================
# MAP OF THE FILE  (in the order run_MCREML.py uses the functions)
# ============================================================================
#   1. read_psam_iids          sample IDs in genotype-row order (the master order)
#   2. read_phenotype          y aligned to that order (NaN = no value)
#   3. build_covariate_table   covariates rebuilt from the raw covariate file
#   4. select_samples          fit samples = have y AND covariates; returns rows, y, X
#   5. read_genotype_chroms    int8 hard calls G (n x m) for the fit samples
#   6. fill_missing_rounded_mean   -9 -> rounded SNP mean, in place (G stays int8)
#   7. genetic_kernels         B -> [K_A B, K_D B]  (one pass over G)
#   8. read_pvar_ids           helper: variant IDs of one .pvar
#   9. read_gene_sets          gene -> list of G columns
#  10. gxg_kernel              B -> W B  (within-gene AxA, low-rank form)
#
#   additive_design / dominance_design are the OLD dense float64 designs.
#   NOTE: the current pipeline (genetic_kernels + gxg_kernel) does not call them;
#   they are kept for reference / comparison only.
#
# A fact used everywhere below: all three kernels are INVARIANT to which allele
# is counted.  Flipping SNP j (g -> 2 - g) flips the sign of z_j, but
#   K_A uses z_j z_j'                     -> unchanged,
#   K_D code (-p/q, 1, -q/p) maps onto itself with p <-> q,
#   GxG uses (z_j o z_k)(z_j o z_k)'      -> unchanged,
#   c-hat uses R_ab s_a s_b, two sign flips -> unchanged.
# So it does not matter whether pgenlib counts REF or ALT alleles.
# ============================================================================
#
# Original header:
#
# Functions that read and align the inputs of MC-REML for one trait:
#   samples (.psam), phenotype (.phen), covariates (rebuilt from the raw
#   covariate files), genotype (.pgen, one or more chromosomes), and the
#   designs built from it (Z for 'A' and 'GxG', Zd for 'D', gene sets for 'GxG').
#
# NOTE: the header line above is partly stale -- Z and Zd are no longer built;
# genetic_kernels / gxg_kernel standardize on the fly from int8 G.
#
# Sample order: the pgen row order (= .psam order) is the order of every
# output, so y[i], X[i, :] and G[i, :] are always the same person.
#
# Paths and parameters are set by the run script; nothing here writes data to disk.

import numpy as np
import pandas as pd
import pgenlib
from scipy.sparse.linalg import svds
from scipy.sparse.linalg import LinearOperator      # NOTE: unused in this file
from concurrent.futures import ThreadPoolExecutor   # threads for the kernels;
                                                    # numpy/BLAS release the GIL,
                                                    # so threads run in parallel


# ============================================================================
# 1. SAMPLES
# ============================================================================
# The .psam file lists the samples of a .pgen, one per row, in the SAME order as
# the genotype rows.  This order is the master order of the whole pipeline:
# phenotype and covariates are re-ordered to match it, never the other way.
def read_psam_iids(pgen_prefix):
    """IIDs in pgen row order, from <pgen_prefix>.psam.

    The .psam header starts with '#' ('#FID IID ...' or '#IID ...'); the '#' is
    stripped so the IID column is found either way.  Raises if an IID appears
    twice, since every other input is matched to the genotype by IID.
    """
    # Tab-separated; the first line is the header, e.g. "#FID  IID  SEX".
    psam = pd.read_csv(pgen_prefix + '.psam', sep='\t')
    # '#FID' -> 'FID' or '#IID' -> 'IID', so psam['IID'] works for both layouts.
    psam.columns = [c.lstrip('#') for c in psam.columns]
    # Plain numpy array of IIDs; position i = genotype row i.
    iids = psam['IID'].to_numpy()
    # Matching by IID only works if IIDs are unique; a duplicate would make
    # reindex() below ambiguous, so fail early.
    if len(set(iids)) != len(iids):
        raise ValueError('duplicate IIDs in %s.psam' % pgen_prefix)
    return iids


# ============================================================================
# 2. PHENOTYPE
# ============================================================================
def read_phenotype(phen_file, trait, iids):
    """Phenotype aligned to iids (pgen order): a float array of len(iids),
    NaN where the person has no value (missing, or absent from the file).

    phen_file is tab-separated with columns FID, IID, <trait> (Phenotype/phen/<trait>.phen).
    """
    phen = pd.read_csv(phen_file, sep='\t')
    # set_index('IID')   : look rows up by IID instead of by position
    # [trait]            : keep the trait column only (a Series indexed by IID)
    # .reindex(iids)     : re-order to the genotype order; an IID that is not in
    #                      the phenotype file gets NaN
    # .to_numpy(float)   : plain float array of length len(iids)
    # NOTE: reindex needs the IIDs to have the same TYPE in both files (e.g. both
    # int).  If one is read as str and the other as int, nothing matches and the
    # check below raises -- which is the intended loud failure.
    y_all = phen.set_index('IID')[trait].reindex(iids).to_numpy(dtype=float)
    # If not a single genotype sample has a value, the files do not belong together.
    if np.isnan(y_all).all():
        raise ValueError('no sample of %s matches the genotype IIDs' % phen_file)
    return y_all


# ============================================================================
# 3. COVARIATES
# ============================================================================
# The covariate matrix X is rebuilt from the raw covariate file every run (not
# stored), with the same rules as get_phen.ipynb.  Columns:
#   quantitative : age, age^2, PC1..PC20 (+ dilution factor for biochemistry)
#   categorical  : assessment center, genotyping array, sex -> one-hot dummies
#   intercept    : 'const' = 1
def build_covariate_table(cov_file, trait):
    """Covariate table for all samples of cov_file, rebuilt from the raw file with
    the rules of get_phen.ipynb (verified against Phenotype/resid_phen).

    Rows: samples (index = IID).  Columns: age, age^2, PC1-20, one-hot assessment
    center / genotyping array / sex (first level dropped), intercept 'const';
    biochemistry traits (name starts with 'biochemistry_') also get
    biochemistry_dilutionfactor.  Samples with a missing quantitative covariate
    are dropped.  Nearly constant columns are removed later, on the fit samples
    (select_samples).
    """
    # Categorical covariates (each becomes several 0/1 dummy columns).
    cov_list = ['cov_ASSESS_CENTER', 'cov_GENO_ARRAY', 'cov_SEX']
    # Quantitative covariates (used as they are, one column each).
    qcov_list = ['cov_AGE', 'cov_AGE_SQ'] + ['PC%d' % x for x in np.arange(1, 21)]
    # Biochemistry traits are also adjusted for the sample dilution factor.
    if trait.startswith('biochemistry_'):
        qcov_list = qcov_list + ['biochemistry_dilutionfactor']

    df_cov = pd.read_csv(cov_file, sep='\t')
    # Guard against pairing a main trait with the biochemistry covariate file
    # (that file has the extra column; its sample set can differ).
    if not trait.startswith('biochemistry_') and 'biochemistry_dilutionfactor' in df_cov.columns:
        raise ValueError('%s is a main trait but %s is the biochemistry covariate file'
                         % (trait, cov_file))

    # --- Categorical covariates -> one-hot, on ALL rows, first level dropped.
    # drop_first=True: the first level of each factor is the reference (absorbed
    # by the intercept); keeping all levels plus 'const' would make X rank deficient.
    temp_df_cov = df_cov[['FID', 'IID'] + cov_list].copy()
    temp_df_cov.index = temp_df_cov['IID']
    temp_df_cov = pd.get_dummies(temp_df_cov[cov_list], columns=cov_list, drop_first=True)
    # NOTE: a missing categorical value becomes all-zero dummies (= reference
    # level), because get_dummies uses dummy_na=False by default.

    # --- Quantitative covariates; drop samples with any missing value.
    # (For biochemistry traits this also drops samples without a dilution factor.)
    df_reg = df_cov[['FID', 'IID'] + qcov_list].copy()
    df_reg.dropna(inplace=True)
    df_reg.index = df_reg['IID']
    # Left join on IID: keeps only the samples that survived dropna, and adds
    # their dummy columns.
    df_reg = df_reg.join(temp_df_cov)
    # The PCs in the raw file have a tiny scale (unit-norm eigenvectors, entries
    # ~ 1/sqrt(N)).  Multiplying by sqrt(N) puts them on a unit-variance-like
    # scale, as get_phen.ipynb does.  Scaling a column of X does not change the
    # fit (projection M = I - X(X'X)^-1 X' is scale invariant); it only helps
    # numerically.
    PC_list = [x for x in df_reg.columns if 'PC' in x]
    df_reg[PC_list] *= np.sqrt(df_reg.shape[0])
    df_reg['const'] = 1                                    # intercept
    df_reg.drop(['FID', 'IID'], axis=1, inplace=True)      # IDs are the index now
    # float32 to keep the table small; select_samples converts X to float64.
    return df_reg.astype(dtype=np.float32)


# ============================================================================
# 4. FIT SAMPLES
# ============================================================================
def select_samples(iids, y_all, df_reg):
    """Samples of the fit and their y and X, in pgen order.

    Samples: all pgen samples with both the phenotype and the covariates.

    Returns (rows, y, X, x_cols):
      rows   0-based pgen rows, strictly increasing (as pgenlib requires)
      y      phenotype of these samples (float64)
      X      covariates of these samples (float64), nearly constant columns
             (sd < 0.01, except 'const') dropped as in get_phen.ipynb
      x_cols names of the columns of X
    """
    # A genotype row i is used iff
    #   y_all[i] is not NaN          (the person has the phenotype), and
    #   iids[i] is in df_reg.index   (the person has all covariates).
    # np.flatnonzero returns the indices in increasing order, which pgenlib's
    # sample_subset requires and which keeps the pgen order.
    rows = np.flatnonzero(~np.isnan(y_all) & pd.Index(iids).isin(df_reg.index))
    y = y_all[rows]

    # Covariate rows of the same people, in the same (pgen) order.
    temp_df_reg = df_reg.loc[iids[rows]]
    # Drop columns that are (nearly) constant IN THE FIT SAMPLES: e.g. a dummy for
    # an assessment center with almost nobody in this subset.  Such a column
    # carries no information and can make X'X singular.  'const' is constant by
    # design and is kept.
    drop_list = [x for x in temp_df_reg if (temp_df_reg[x].std() < 1e-2) & (x != 'const')]
    temp_df_reg = temp_df_reg.drop(columns=drop_list)
    X = temp_df_reg.to_numpy(dtype=np.float64)
    # X must have full column rank: the REML projection uses (X'X)^-1.
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError('covariate matrix X is rank deficient in the fit samples')
    return rows, y, X, list(temp_df_reg.columns)


# ============================================================================
# 5. GENOTYPE
# ============================================================================
def read_genotype_chroms(pgen_prefixes, rows, iids):
    """Hard calls of all variants of one or more chromosomes, side by side.

    pgen_prefixes: one prefix per chromosome, in the column order wanted.
    Every .psam must list exactly iids (same samples, same order), otherwise
    rows would point to different people in different files.

    Returns (G, snp_chrom, snp_idx):
      G          (len(rows), total m) int8: 0/1/2, -9 missing; row i = pgen row rows[i]
      snp_chrom  for each column, the position of its chromosome in pgen_prefixes
      snp_idx    for each column, its 0-based row in that chromosome's .pvar
    """
    # --- Pass 1: check every chromosome's .psam and count its variants, so the
    # whole buffer can be allocated once (no growing / concatenating copies).
    m_list = []
    for prefix in pgen_prefixes:
        # rows are indices into iids; they only mean the same people in every
        # file if every .psam has exactly the same IIDs in the same order.
        if not np.array_equal(read_psam_iids(prefix), iids):
            raise ValueError('%s.psam does not match iids' % prefix)
        reader = pgenlib.PgenReader((prefix + '.pgen').encode())
        m_list.append(reader.get_variant_ct())             # number of SNPs on this chr
        reader.close()

    # --- Pass 2: read. pgenlib writes variant-major (one row per SNP), so the
    # buffer is (total m, n).  int8 = 1 byte per call: n * m bytes in total
    # (100K x 1.17M ~ 117 GB genome-wide).
    buf = np.empty((sum(m_list), len(rows)), dtype=np.int8)
    start = 0
    for prefix, m in zip(pgen_prefixes, m_list):
        # sample_subset: only the fit samples (rows) are decoded, in that order.
        reader = pgenlib.PgenReader((prefix + '.pgen').encode(),
                                    sample_subset=rows.astype(np.uint32))
        # All variants 0..m-1 of this chromosome into its slice of buf
        # (buf[start:start+m] is a view, so pgenlib writes straight into buf).
        # Values: 0/1/2 allele counts, -9 = missing call.
        reader.read_range(0, m, buf[start:start + m])
        reader.close()
        start += m

    # Column bookkeeping, used later to map gene tables onto G:
    #   snp_chrom[j] = which chromosome (position in pgen_prefixes) column j is on
    #   snp_idx[j]   = row of column j within that chromosome's .pvar
    # e.g. with chr21 (m21 SNPs) + chr22: column m21 has snp_chrom 1, snp_idx 0.
    snp_chrom = np.repeat(np.arange(len(pgen_prefixes)), m_list)
    snp_idx = np.concatenate([np.arange(m) for m in m_list])
    # buf.T is a VIEW (no copy) of shape (n, m).  Because buf is C-ordered
    # (m, n), G = buf.T is Fortran-ordered: each COLUMN (one SNP) is contiguous
    # in memory, which makes the column blocks G[:, a:a+block] used by the
    # kernels cheap to read.
    return buf.T, snp_chrom, snp_idx


# ============================================================================
# 8-9. GENE SETS
# ============================================================================
def read_pvar_ids(pgen_prefix):
    """Variant IDs in .pvar row order (skips '##' meta lines)."""
    # A .pvar may start with '##...' meta lines, then one header line
    # '#CHROM POS ID REF ALT ...', then one line per variant.  Count the '##'
    # lines so read_csv can skip them and use the '#CHROM ...' line as header.
    with open(pgen_prefix + '.pvar') as f:
        n_meta = 0
        for line in f:
            if not line.startswith('##'):
                break                                # reached the '#CHROM' header
            n_meta += 1
    # usecols=['ID'] reads only the ID column (fast, little memory).
    pvar = pd.read_csv(pgen_prefix + '.pvar', sep='\t', skiprows=n_meta, usecols=['ID'])
    return pvar['ID'].to_numpy()                     # position j = .pvar row j


def read_gene_sets(gene_table_files, pgen_prefixes, snp_chrom_kept, snp_idx_kept, min_snps=2):
    """Gene sets for the within-gene AxA kernel, as column indices of Z.  ('GxG')

    gene_table_files: one long-format gene table per chromosome, in the same
    order as pgen_prefixes (columns gene_id, snp_idx, snp_id, ...; snp_idx =
    0-based .pvar row, snp_id = ID at that row).
    snp_chrom_kept / snp_idx_kept: chromosome position and .pvar row of each
    column of Z, i.e. snp_chrom[keep], snp_idx[keep].

    SNPs removed by the MAF filter are dropped from their genes; genes left with
    fewer than min_snps SNPs are dropped (no within-gene pair).

    Returns (gene_ids, gene_sets): gene_sets[g] = Z columns of gene gene_ids[g].
    """
    # NOTE: the docstring still describes the old dense-Z pipeline.  Today it is
    # called with ALL columns, read_gene_sets(tables, prefixes, snp_chrom, snp_idx),
    # so "column of Z" = "column of G" and no SNP is removed here.
    #
    # Gene region definition (decided upstream, in the gene-table script):
    # exons + TSS +-1 kb of every gene.  Introns are NOT included, which is why
    # a gene's columns are usually not one contiguous range.
    #
    # Output, e.g.:  gene_ids[0]  = 'ENSG00000141956.9'
    #                gene_sets[0] = array([312, 313, 318, ...])  -> G[:, gene_sets[0]]
    # Column numbers are GLOBAL (all chromosomes read into G, side by side), so
    # they depend on which chromosomes were read: always build gene_sets from the
    # SAME read as G, never reuse them with another G.
    gene_ids, gene_sets = [], []
    # k = position of the chromosome in pgen_prefixes (0, 1, ...), matching snp_chrom.
    for k, (gt_file, prefix) in enumerate(zip(gene_table_files, pgen_prefixes)):
        # Long format: one row per (gene, SNP); a SNP in two genes has two rows.
        gt = pd.read_csv(gt_file, sep='\t')

        # Safety check: the gene table must have been built on THIS .pvar --
        # the ID at each listed .pvar row must equal the listed snp_id.
        # NOTE: an snp_idx past the end of the .pvar raises IndexError here
        # (still a loud failure, just a different exception type).
        pvar_ids = read_pvar_ids(prefix)
        if not np.array_equal(pvar_ids[gt['snp_idx'].to_numpy()], gt['snp_id'].to_numpy()):
            raise ValueError('%s does not match %s.pvar' % (gt_file, prefix))

        # Lookup table for this chromosome: .pvar row -> global column of G
        # (-1 = that .pvar row is not a column of G).
        col_of = np.full(len(pvar_ids), -1, dtype=np.int64)
        on_k = np.flatnonzero(snp_chrom_kept == k)        # global columns on chr k
        col_of[snp_idx_kept[on_k]] = on_k                 # their .pvar rows -> columns
        # e.g. chr22 read after chr21 (m21 SNPs): .pvar row 0 of chr22 -> column m21.

        # Translate every gene-table row to its G column; drop rows not in G.
        gt['col'] = col_of[gt['snp_idx'].to_numpy()]
        gt = gt[gt['col'] >= 0]
        # One entry per gene.  sort=False keeps the order of first appearance in
        # the table (the tables are ordered by Ensembl ID, not by position --
        # irrelevant for the kernel, which sums over genes).
        for gene_id, df_g in gt.groupby('gene_id', sort=False):
            # A gene needs >= 2 SNPs to have at least one pair (j < k).
            if len(df_g) >= min_snps:
                gene_ids.append(gene_id)
                gene_sets.append(df_g['col'].to_numpy())
    return gene_ids, gene_sets


# ============================================================================
# 10. GxG KERNEL  (within-gene additive x additive, low-rank form)
# ============================================================================
# Target matrix (exact, r = full rank):
#     W = 1/(P c) * sum_g sum_{j<k in g} h_jk h_jk',   h_jk = z_j o z_k
# Identity used:  sum_{j<k} h_jk h_jk' = 1/2 [ K_g o K_g - D_g D_g' ],
#     K_g = Z_g Z_g',  D_g = Z_g o Z_g  (the j = k terms of K_g o K_g).
# Low rank: replace ONE factor K_g by its rank-r truncation Kt_g = sum_s lam_s q_s q_s'.
# Factorization:  Kt_g o K_g = A_g A_g',  A_g = [ sqrt(lam_s) q_s o z_j ]_{s, j}
# because (q q') o (z z') = (q o z)(q o z)'.
# A_g has r_g * m_g columns and is built per gene on every apply (never stored).
def gxg_kernel(G, gene_sets, r=9999, c=None, n_threads=1, max_cols=512):
    """Within-gene additive-by-additive kernel as a function B -> W B.  ('GxG')

        W = 1/(2 P c) * sum_g [ Kt_g o K_g - D_g D_g' ]
        K_g  = Z_g Z_g',  Kt_g = rank-r truncation of K_g,  D_g = Z_g o Z_g
    With r >= rank(Z_g), Kt_g = K_g and gene g gives exactly 2 * sum_{j<k} h_jk h_jk',
    h_jk = z_j o z_k (the j = k terms are removed by D_g D_g').
    Z_g = HWE-standardized G[:, gene_sets[g]], (g - 2p) / sqrt(2pq), rebuilt from G
    on every apply; nothing n-sized is stored.  With eigenpairs (lam_s, v_s) of
    Z_g' Z_g and q_s = Z_g v_s / sqrt(lam_s), Kt_g = sum_s lam_s q_s q_s'.
    Per gene, A_g = [sqrt(lam_s) q_s o Z_g]_s is built on every apply and
        (Kt_g o K_g) B = A_g (A_g' B)
    P = sum_g m_g (m_g - 1) / 2 over the used SNPs.
    c: divisor; None = c-hat = 1 + (1/P) sum_g sum_{a<b} R_ab s_a s_b,
       R_g = Z_g' Z_g / n, s_a = mean(z_a^3)  (third-moment plug-in for mean Var(z_a z_b)).
    Monomorphic SNPs are dropped from their genes; genes left with < 2 SNPs are dropped.

    G: (n, m) int8 with no missing calls.  gene_sets: from read_gene_sets on this G.
    max_cols caps the width of each temporary A block, (n, s * m_g).
    Returns K: K(B) gives (n,) or (n, k).  K.P, K.n_genes, K.rank (r_g per gene), K.c.
    """
    n = G.shape[0]

    # ------------------------------------------------------------------ setup
    # Runs ONCE per gene.  Keeps only small per-gene objects (size m_g or
    # m_g x r_g), never anything of length n.
    def setup(cols):
        cols = np.asarray(cols)
        p = G[:, cols].mean(axis=0, dtype=np.float64) / 2
        ok = (p > 0) & (p < 1)                               # drop monomorphic SNPs
        cols, p = cols[ok], p[ok]
        if len(cols) < 2:
            return None
        ctr = (2 * p).astype(np.float32)
        scl = (1 / np.sqrt(2 * p * (1 - p))).astype(np.float32)
        Z = (G[:, cols].astype(np.float32) - ctr) * scl
        Zd = Z.astype(np.float64)
        mg = Zd.shape[1]

        # c-hat ingredient in O(n m_g), Z_g' Z_g never formed:
        #   sum_{a != b} R_ab s_a s_b = s' R s - sum_a R_aa s_a^2
        #                             = (||Z s||^2 - sum_a ||z_a||^2 s_a^2) / n
        s = (Zd * Zd * Zd).mean(axis=0)                      # per-SNP skewness
        Zs = Zd @ s                                          # one mat-vec
        cross = (Zs @ Zs - (Zd * Zd).sum(axis=0) @ (s * s)) / n

        # Leading r singular triplets of Z_g: lam_s = sigma_s^2 are the eigenvalues
        # of K_g = Z Z', v_s the right singular vectors.  Truncated SVD (ARPACK)
        # costs O(n m_g r) per iteration; ARPACK needs k < min(n, m_g), so when r
        # reaches full rank the gene is factorized by an exact thin SVD instead.
        if r >= mg - 1:
            _, sv, Vt = np.linalg.svd(Zd, full_matrices=False)  # descending
        else:
            _, sv, Vt = svds(Zd, k=r, v0=np.ones(mg))           # fixed v0: reproducible
            idx = np.argsort(-sv)                                # svds: no fixed order
            sv, Vt = sv[idx], Vt[idx]
        lam, V = sv ** 2, Vt.T
        rg = min(r, int((lam > 1e-8 * lam[0]).sum()))        # skip null directions
        lam, V = lam[:rg], V[:, :rg]
        T = (V / np.sqrt(lam)).astype(np.float32)            # q = Z T has unit columns
        return cols, ctr, scl, T, lam.astype(np.float32), cross

    # Setup of all genes in parallel; genes that returned None are left out.
    with ThreadPoolExecutor(n_threads) as ex:
        genes = [g for g in ex.map(setup, gene_sets) if g is not None]
    if not genes:
        raise ValueError("no gene with >= 2 polymorphic SNPs")
    m_g = np.array([len(g[0]) for g in genes])               # SNPs per gene (after filter)
    r_g = np.array([g[3].shape[1] for g in genes])           # kept rank per gene
    # Total number of within-gene pairs j < k: the averaging constant of W
    # (plays the role m plays in K_A = Z Z' / m).
    P = int((m_g * (m_g - 1) // 2).sum())
    if c is None:
        
        c = 1.0 + sum(g[5] for g in genes) / (2.0 * P)
        # c-hat is an estimate and could in principle be <= 0; dividing by it
        # would flip the sign of the GxG component, so refuse.
        if not c > 0:
            raise ValueError(f"c-hat = {c!r} is not positive; check the SNP filter")

    # Load balancing: cost of a gene per apply ~ m_g * (r_g + 1) columns
    # (A_g: m_g r_g, D_g: m_g).  Sort genes from most to least expensive and deal
    # them round-robin to the threads, so every thread gets a similar total.
    order = np.argsort(-(m_g * (r_g + 1)))
    chunks = [order[t::n_threads] for t in range(n_threads)]

    # ------------------------------------------------------------------ apply
    # K(B) = W B for B of shape (n,) or (n, k); called every CG iteration.
    def K(B):
        B2 = np.asarray(B, dtype=np.float64).reshape(n, -1)   # (n,) -> (n, 1)
        k = B2.shape[1]
        B32 = B2.astype(np.float32)                           # float32 for the GEMMs

        def part(idx):
            # Each thread accumulates the contribution of its own genes.
            acc = np.zeros((n, k))                            # float64 accumulator
            for i in idx:
                cols, ctr, scl, T, lam, _ = genes[i]
                # Rebuild Z_g from int8 G (same formula as in setup).
                Z = (G[:, cols].astype(np.float32) - ctr) * scl          # (n, m_g)
                mg = Z.shape[1]
                # Qs[:, s] = sqrt(lam_s) q_s = Z (v_s / sqrt(lam_s)) sqrt(lam_s) = Z v_s.
                # (So Qs equals Z @ V; going through T keeps one stored factor.)
                Qs = (Z @ T) * np.sqrt(lam)                              # (n, r_g)
                # Number of components s per A block so that the block has at most
                # max_cols = 512 columns: n x 512 float32 ~ 200 MB at n = 100K.
                s_step = max(1, max_cols // mg)
                for a in range(0, Qs.shape[1], s_step):
                    Qa = Qs[:, a:a + s_step]                             # (n, ns)
                    ns = Qa.shape[1]
                    # Broadcast (n, ns, 1) * (n, 1, m_g) -> (n, ns, m_g), then
                    # flatten to (n, ns*m_g).  Column (s, j) = sqrt(lam_s) q_s o z_j.
                    # A A' = sum_s lam_s (q_s q_s') o (sum_j z_j z_j') = (part of Kt_g) o K_g.
                    A = (Qa[:, :, None] * Z[:, None, :]).reshape(n, ns * mg)
                    # Two GEMMs: (A' B) is (ns*m_g, k), then A (.) is (n, k).
                    acc += A @ (A.T @ B32)
                # Remove the j = k terms of the Hadamard square:
                # D_g D_g' = sum_j (z_j o z_j)(z_j o z_j)'.
                # (With truncation, the A part is approximate but this term is exact,
                # as in the original stored-A code.)
                D = Z * Z
                acc -= D @ (D.T @ B32)
            return acc

        # Sum the threads' partial results, then the scaling 1 / (2 P c):
        # the 2 undoes A A' - D D' = 2 * sum_{j<k} h h'.
        with ThreadPoolExecutor(n_threads) as ex:
            out = sum(ex.map(part, chunks)) / (2.0 * P * c)
        # Same shape as the input: (n,) in -> (n,) out.
        return out.ravel() if np.ndim(B) == 1 else out

    # K is a closure: it "remembers" genes, P, c, chunks from the setup above.
    # Diagnostics attached as attributes (written to the result table):
    K.P, K.n_genes, K.rank, K.c = P, len(genes), r_g, c
    return K


# ============================================================================
# 6. MISSING CALLS
# ============================================================================
def fill_missing_rounded_mean(G, block=4096):
    """In place: each missing call (-9) of G becomes the rounded mean of its SNP.

    int8 cannot hold the exact mean; the rounded mean keeps G integer.  Done
    block by block so the temporary arrays stay small.
    """
    # Must run BEFORE genetic_kernels / gxg_kernel: those assume G has only 0/1/2.
    # (A -9 left in G would be used as a genotype value; in the dominance lookup
    # tab[:, j][-9] it would even raise IndexError.)
    n, m = G.shape
    for a in range(0, m, block):
        g = G[:, a:a + block]                            # VIEW into G: writes go to G
        obs = g >= 0
        # Mean genotype of the observed calls, per SNP (int64 sum avoids overflow).
        mu0 = np.where(obs, g, 0).sum(axis=0, dtype=np.int64) / obs.sum(axis=0)
        # NOTE: a SNP with ALL calls missing gives 0/0 = NaN here; rint(NaN) cast
        # to int8 is undefined.  Not an issue after upstream QC, but worth knowing.
        i_na, j_na = np.nonzero(~obs)                    # positions of missing calls
        # Rounded mean: e.g. mean 0.6 -> 1.  For common SNPs this is usually the
        # heterozygote, so heavy missingness slightly inflates heterozygosity
        # (a small departure from HWE; one source of the c-hat vs exact gap).
        g[i_na, j_na] = np.rint(mu0[j_na]).astype(np.int8)
    return G


# ============================================================================
# 7. ADDITIVE + DOMINANCE KERNELS  (one pass over G)
# ============================================================================
# Math (additive), with Omega = diag(w), w_j = 1/sqrt(2 p_j q_j), mu_j = 2 p_j:
#     Z = (G - 1 mu') Omega,     K_A = Z Z' / m_used.
# Z is never formed.  Per SNP block b:
#     C_b   = Z_b' B = Omega_b (G_b' B - mu_b (1' B))
#     Z_b C_b = G_b (Omega_b C_b) - 1 (mu_b' Omega_b C_b)
# so with WC_b = Omega_b^2 (G_b' B - mu_b (1' B)):
#     Z_b Z_b' B = G_b WC_b - 1 (mu_b' WC_b)            (rank-one correction)
# and K_A B = sum_b [ G_b WC_b - 1 (mu_b' WC_b) ] / m_used.
#
# Math (dominance), GCTA / Zhu et al. code with p = mu/2, q = 1 - p:
#     g = 0 -> -p/q,  g = 1 -> 1,  g = 2 -> -q/p   (mean 0, var 1 under HWE)
#     K_D = D D' / m_used,   K_D B = sum_b D_b (D_b' B) / m_used.
def genetic_kernels(G, components=('A', 'D'), block=512, n_threads=1):
    """Additive and dominance kernels in ONE pass over G: B -> [K_c B for c in components].

    'A': K_A = Z Z' / m,  Z = (G - 2p) / sqrt(2pq)   (HWE standardization, p = mean / 2)
    'D': K_D = D D' / m,  D = GCTA dominance code with p = mu / 2, q = 1 - p:
             g = 0 -> -p/q,   g = 1 -> 1,   g = 2 -> -q/p
         (mean 0 and variance 1 under HWE, so no further standardization).
    Per SNP block, G_b is cast to float32 once and shared by both kernels; D_b is
    built by a per-SNP table lookup.  Blocks are split over n_threads threads (set
    BLAS to 1 thread).  Monomorphic SNPs (p = 0 or 1) get weight 0 in both kernels.

    G: (n, m) int8 with no missing calls (fill_missing_rounded_mean first).
    Returns K: K(B) gives a list of (n,) or (n, k) arrays, one per component.
    """
    # Only 'A' and 'D' are handled here; GxG is gxg_kernel (separate, gene blocks).
    if not set(components) <= {'A', 'D'} or len(components) == 0:
        raise ValueError(f"components must be a non-empty subset of ('A', 'D'), got {components}")
    n, m = G.shape
    # Block start columns 0, 512, 1024, ...; thread t gets blocks t, t+T, t+2T, ...
    # (round-robin, so all threads get about the same number of blocks).
    starts = list(range(0, m, block))
    chunks = [starts[t::n_threads] for t in range(n_threads)]   # blocks of each thread

    # --- SNP means (one pass, threaded).  Each thread writes its own slices of
    # mu, so there is no write conflict.
    mu = np.empty(m)

    def stats(blocks):
        for a in blocks:
            mu[a:a + block] = G[:, a:a + block].mean(axis=0, dtype=np.float64)

    with ThreadPoolExecutor(n_threads) as ex:
        list(ex.map(stats, chunks))                      # list() waits for all threads

    # --- Per-SNP constants (HWE scaling).
    used = (mu > 0) & (mu < 2)                           # not monomorphic: 0 < p < 1
    m_used = int(used.sum())                             # divisor of both kernels
    p = np.where(used, mu / 2, 0.5)                      # 0.5 only avoids 0/0 below
    q = 1 - p
    # Additive weight w_j = 1/sqrt(2 p q); 0 for monomorphic SNPs, which makes
    # their column of Z zero, i.e. they drop out of K_A.
    w = np.zeros(m)                                      # additive: 1 / sqrt(2pq)
    w[used] = 1.0 / np.sqrt(2 * p[used] * q[used])
    # Dominance lookup table, one column per SNP: row g holds the code of genotype g.
    # Monomorphic SNPs keep all-zero codes -> they drop out of K_D.
    tab = np.zeros((3, m), dtype=np.float32)             # dominance code for g = 0, 1, 2
    tab[0, used] = -p[used] / q[used]
    tab[1, used] = 1.0
    tab[2, used] = -q[used] / p[used]

    do_A = 'A' in components
    do_D = 'D' in components

    # ------------------------------------------------------------------ apply
    def K(B):
        B2 = np.asarray(B, dtype=np.float64).reshape(n, -1)   # (n,) -> (n, 1)
        k = B2.shape[1]
        B32 = B2.astype(np.float32)                           # float32 for the GEMMs
        colsum = B2.sum(axis=0)                          # 1' B  (needed by the centering)

        def part(blocks):
            # Per-thread accumulators (float64 to limit rounding error over ~2K blocks).
            acc_A = np.zeros((n, k)) if do_A else None   # sum_b G_b WC_b
            corr_A = np.zeros(k) if do_A else None       # sum_b mu_b' WC_b
            acc_D = np.zeros((n, k)) if do_D else None   # sum_b D_b (D_b' B)
            # Reusable buffer for one block of dominance codes (Fortran order so
            # that writing one column d[:, j] is contiguous).
            d = np.empty((n, block), dtype=np.float32, order='F') if do_D else None
            for a in blocks:
                Gb = G[:, a:a + block]                   # view, (n, nb) int8
                nb = Gb.shape[1]                         # last block may be < 512
                if do_A:
                    g = Gb.astype(np.float32)            # the ONE cast of this block
                    wb, mb = w[a:a + block], mu[a:a + block]
                    # WC = Omega^2 (G_b' B - mu_b (1' B)),  shape (nb, k)
                    #   g.T @ B32                     : G_b' B
                    #   mb[:, None] * colsum[None, :] : mu_b (1' B)   (centering)
                    #   (wb * wb)[:, None]            : Omega_b^2     (scaling, twice)
                    WC = (g.T @ B32 - mb[:, None] * colsum[None, :]) * (wb * wb)[:, None]
                    acc_A += g @ WC.astype(np.float32)   # G_b WC_b   (n, k)
                    corr_A += mb @ WC                    # mu_b' WC_b (k,)
                if do_D:
                    # Build D_b by lookup: for SNP j, tab[:, a+j] is its 3 codes and
                    # Gb[:, j] (values 0/1/2) picks one code per person.
                    for j in range(nb):                  # per-SNP lookup of the code
                        d[:, j] = tab[:, a + j][Gb[:, j]]
                    db = d[:, :nb]
                    acc_D += db @ (db.T @ B32)           # D_b (D_b' B)
            return acc_A, corr_A, acc_D

        with ThreadPoolExecutor(n_threads) as ex:
            parts = list(ex.map(part, chunks))           # one (acc_A, corr_A, acc_D) per thread
        out = []
        # Results in the order of `components`, e.g. ['A', 'D'] -> [K_A B, K_D B].
        for c in components:
            if c == 'A':
                # sum over threads, then the rank-one correction:
                # K_A B = (sum G_b WC_b - 1 (sum mu_b' WC_b)) / m_used
                r = (sum(t[0] for t in parts) - sum(t[1] for t in parts)[None, :]) / m_used
            else:
                r = sum(t[2] for t in parts) / m_used
            out.append(r.ravel() if np.ndim(B) == 1 else r)
        return out

    return K