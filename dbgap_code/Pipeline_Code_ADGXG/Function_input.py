import numpy as np
import pandas as pd
import pgenlib
from scipy.sparse.linalg import svds
from concurrent.futures import ThreadPoolExecutor
import scipy.sparse.linalg as spla

def read_psam_iids(pgen_prefix):
    """IIDs in pgen row order, from <pgen_prefix>.psam.

    The .psam header starts with '#' ('#FID IID ...' or '#IID ...'); the '#' is
    stripped so the IID column is found either way.  Raises if an IID appears
    twice, since every other input is matched to the genotype by IID.
    """
    psam = pd.read_csv(pgen_prefix + '.psam', sep='\t')
    psam.columns = [c.lstrip('#') for c in psam.columns]
    iids = psam['IID'].to_numpy()
    if len(set(iids)) != len(iids):
        raise ValueError('duplicate IIDs in %s.psam' % pgen_prefix)
    return iids


def read_phenotype(phen_file, trait, iids):
    """Phenotype aligned to iids (pgen order): a float array of len(iids),
    NaN where the person has no value (missing, or absent from the file).

    phen_file is tab-separated with columns FID, IID, <trait> (Phenotype/phen/<trait>.phen).
    """
    phen = pd.read_csv(phen_file, sep='\t')
    y_all = phen.set_index('IID')[trait].reindex(iids).to_numpy(dtype=float)
    if np.isnan(y_all).all():
        raise ValueError('no sample of %s matches the genotype IIDs' % phen_file)
    return y_all


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
    cov_list = ['cov_ASSESS_CENTER', 'cov_GENO_ARRAY', 'cov_SEX']
    qcov_list = ['cov_AGE', 'cov_AGE_SQ'] + ['PC%d' % x for x in np.arange(1, 21)]
    if trait.startswith('biochemistry_'):
        qcov_list = qcov_list + ['biochemistry_dilutionfactor']

    df_cov = pd.read_csv(cov_file, sep='\t')
    if not trait.startswith('biochemistry_') and 'biochemistry_dilutionfactor' in df_cov.columns:
        raise ValueError('%s is a main trait but %s is the biochemistry covariate file'
                         % (trait, cov_file))

    temp_df_cov = df_cov[['FID', 'IID'] + cov_list].copy()
    temp_df_cov.index = temp_df_cov['IID']
    temp_df_cov = pd.get_dummies(temp_df_cov[cov_list], columns=cov_list, drop_first=True)

    df_reg = df_cov[['FID', 'IID'] + qcov_list].copy()
    df_reg.dropna(inplace=True)
    df_reg.index = df_reg['IID']
    df_reg = df_reg.join(temp_df_cov)
    PC_list = [x for x in df_reg.columns if 'PC' in x]
    df_reg[PC_list] *= np.sqrt(df_reg.shape[0])
    df_reg['const'] = 1
    df_reg.drop(['FID', 'IID'], axis=1, inplace=True)
    return df_reg.astype(dtype=np.float32)


def select_samples(iids, y_all, df_reg):
    """Samples of the fit and their y and X, in pgen order.

    Samples: all pgen samples with both the phenotype and the covariates.

    Returns (rows, y, X, x_cols):
      rows   0-based pgen rows, strictly increasing (as pgenlib requires)
      y      phenotype of these samples (float64)
      X      covariates of these samples (float64): nearly constant columns
             (sd < 0.01, except 'const') dropped as in get_phen.ipynb, then
             columns that are linear combinations of others dropped (e.g. the
             dummies of a factor whose reference level has nobody in the fit
             samples sum to 'const').  The fit only uses the column space of X,
             so dropping them does not change it.
      x_cols names of the columns of X
    """
    rows = np.flatnonzero(~np.isnan(y_all) & pd.Index(iids).isin(df_reg.index))
    y = y_all[rows]

    temp_df_reg = df_reg.loc[iids[rows]]
    drop_list = [x for x in temp_df_reg if (temp_df_reg[x].std() < 1e-2) & (x != 'const')]
    temp_df_reg = temp_df_reg.drop(columns=drop_list)

    # Drop linearly dependent columns: Gram-Schmidt over 'const' first, then the
    # other columns in table order; a column is kept only if it adds a new direction.
    order = ['const'] + [c for c in temp_df_reg.columns if c != 'const']
    Q = np.empty((len(rows), 0))
    keep = set()
    for c in order:
        v = temp_df_reg[c].to_numpy(dtype=np.float64)
        r = v - Q @ (Q.T @ v)
        r = r - Q @ (Q.T @ r)                        # second pass: keeps Q orthonormal
        if np.linalg.norm(r) > 1e-8 * np.linalg.norm(v):
            Q = np.column_stack([Q, r / np.linalg.norm(r)])
            keep.add(c)
    temp_df_reg = temp_df_reg[[c for c in temp_df_reg.columns if c in keep]]

    X = temp_df_reg.to_numpy(dtype=np.float64)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError('covariate matrix X is rank deficient in the fit samples')
    return rows, y, X, list(temp_df_reg.columns)

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
    m_list = []
    for prefix in pgen_prefixes:
        if not np.array_equal(read_psam_iids(prefix), iids):
            raise ValueError('%s.psam does not match iids' % prefix)
        reader = pgenlib.PgenReader((prefix + '.pgen').encode())
        m_list.append(reader.get_variant_ct())
        reader.close()

    buf = np.empty((sum(m_list), len(rows)), dtype=np.int8)
    start = 0
    for prefix, m in zip(pgen_prefixes, m_list):
        reader = pgenlib.PgenReader((prefix + '.pgen').encode(),
                                    sample_subset=rows.astype(np.uint32))
        reader.read_range(0, m, buf[start:start + m])
        reader.close()
        start += m

    snp_chrom = np.repeat(np.arange(len(pgen_prefixes)), m_list)
    snp_idx = np.concatenate([np.arange(m) for m in m_list])
    return buf.T, snp_chrom, snp_idx


def read_pvar_ids(pgen_prefix):
    """Variant IDs in .pvar row order (skips '##' meta lines)."""
    with open(pgen_prefix + '.pvar') as f:
        n_meta = 0
        for line in f:
            if not line.startswith('##'):
                break
            n_meta += 1
    pvar = pd.read_csv(pgen_prefix + '.pvar', sep='\t', skiprows=n_meta, usecols=['ID'])
    return pvar['ID'].to_numpy()


def read_gene_sets(gene_table_files, pgen_prefixes, snp_chrom_kept, snp_idx_kept, min_snps=2):
    """Gene sets for the within-gene AxA kernel, as column indices of G.  ('GxG')

    gene_table_files: one long-format gene table per chromosome, in the same
    order as pgen_prefixes (columns gene_id, snp_idx, snp_id, ...; snp_idx =
    0-based .pvar row, snp_id = ID at that row).
    snp_chrom_kept / snp_idx_kept: chromosome position and .pvar row of each
    column of G, i.e. snp_chrom and snp_idx from read_genotype_chroms on the
    same read as G.

    Gene-table SNPs that are not columns of G are dropped from their genes;
    genes left with fewer than min_snps SNPs are dropped (no within-gene pair).

    Returns (gene_ids, gene_sets): gene_sets[g] = G columns of gene gene_ids[g].
    """
    gene_ids, gene_sets = [], []
    for k, (gt_file, prefix) in enumerate(zip(gene_table_files, pgen_prefixes)):
        gt = pd.read_csv(gt_file, sep='\t')

        pvar_ids = read_pvar_ids(prefix)
        if not np.array_equal(pvar_ids[gt['snp_idx'].to_numpy()], gt['snp_id'].to_numpy()):
            raise ValueError('%s does not match %s.pvar' % (gt_file, prefix))

        col_of = np.full(len(pvar_ids), -1, dtype=np.int64)
        on_k = np.flatnonzero(snp_chrom_kept == k)
        col_of[snp_idx_kept[on_k]] = on_k

        gt['col'] = col_of[gt['snp_idx'].to_numpy()]
        gt = gt[gt['col'] >= 0]
        for gene_id, df_g in gt.groupby('gene_id', sort=False):
            if len(df_g) >= min_snps:
                gene_ids.append(gene_id)
                gene_sets.append(df_g['col'].to_numpy())
    return gene_ids, gene_sets


def gxg_kernel(G, gene_sets, r=9999, c=None, n_threads=1, max_cols=512):
    """Within-gene additive-by-additive kernel as a LinearOperator W.  ('GxG')

        W = 1/(2 P c) * sum_g [ Kt_g o K_g - Dt_g Dt_g' ]
        K_g = Z_g Z_g',  Zt_g = rank-r truncation of Z_g (truncated SVD),  Kt_g = Zt_g Zt_g',
        Dt_g = Zt_g o Z_g
    Kt_g o K_g = sum_{j,k} g_jk g_jk' with g_jk = zt_j o z_k, and Dt_g Dt_g' is exactly its
    j = k part, so each gene gives sum_{j != k} g_jk g_jk' >= 0: W is PSD for every r.
    With r >= m_g, Zt_g = Z_g, Dt_g = Z_g o Z_g, and gene g gives exactly
    2 * sum_{j<k} h_jk h_jk', h_jk = z_j o z_k.
    Z_g = HWE-standardized G[:, gene_sets[g]], (g - 2p) / sqrt(2pq), rebuilt from G
    on every apply; nothing n-sized is stored.  With the leading r singular
    triplets of Z_g (Z_g' Z_g is never formed), lam_s = sigma_s^2, right singular
    vectors v_s, q_s = Z_g v_s / sqrt(lam_s): Kt_g = sum_s lam_s q_s q_s' and
    Zt_g = sum_s sqrt(lam_s) q_s v_s' = (Z_g V_r) V_r'.
    Per gene, A_g = [sqrt(lam_s) q_s o Z_g]_s is built on every apply and
        (Kt_g o K_g) B = A_g (A_g' B)
    P = sum_g m_g (m_g - 1) / 2 over the used SNPs.
    c: divisor; None = c-hat = 1 + (1/P) sum_g sum_{a<b} R_ab s_a s_b,
       R_g = Z_g' Z_g / n, s_a = mean(z_a^3)  (third-moment plug-in for mean Var(z_a z_b)).
    Monomorphic SNPs are dropped from their genes; genes left with < 2 SNPs are dropped.

    G: (n, m) int8 with no missing calls.  gene_sets: from read_gene_sets on this G.
    max_cols caps the width of each temporary A block, (n, s * m_g).
    Returns W: (n, n) LinearOperator, symmetric (rmatvec = matvec); W @ B for (n,) or (n, k).
    Attributes W.P, W.n_genes, W.rank (r_g per gene), W.c.
    """
    n = G.shape[0]

    def setup(cols):
        cols = np.asarray(cols)
        p = G[:, cols].mean(axis=0, dtype=np.float64) / 2
        ok = (p > 0) & (p < 1)
        cols, p = cols[ok], p[ok]
        if len(cols) < 2:
            return None
        ctr = (2 * p).astype(np.float32)
        scl = (1 / np.sqrt(2 * p * (1 - p))).astype(np.float32)
        Z = (G[:, cols].astype(np.float32) - ctr) * scl
        Zd = Z.astype(np.float64)
        mg = Zd.shape[1]

        s = (Zd * Zd * Zd).mean(axis=0)
        Zs = Zd @ s
        cross = (Zs @ Zs - (Zd * Zd).sum(axis=0) @ (s * s)) / n

        if r >= mg - 1:                                 # svds needs k < m_g: use the full SVD
            _, sv, Vt = np.linalg.svd(Zd, full_matrices=False)
        else:
            _, sv, Vt = svds(Zd, k=r, v0=np.ones(mg))
            idx = np.argsort(-sv)
            sv, Vt = sv[idx], Vt[idx]
        lam, V = sv ** 2, Vt.T
        rg = min(r, int((lam > 1e-8 * lam[0]).sum()))
        lam, V = lam[:rg], V[:, :rg]
        T = (V / np.sqrt(lam)).astype(np.float32)
        full = r >= mg                                  # no truncation: Zt_g = Z_g
        return cols, ctr, scl, T, lam.astype(np.float32), cross, full

    with ThreadPoolExecutor(n_threads) as ex:
        genes = [g for g in ex.map(setup, gene_sets) if g is not None]
    if not genes:
        raise ValueError("no gene with >= 2 polymorphic SNPs")
    m_g = np.array([len(g[0]) for g in genes])
    r_g = np.array([g[3].shape[1] for g in genes])
    P = int((m_g * (m_g - 1) // 2).sum())
    if c is None:
        c = 1.0 + sum(g[5] for g in genes) / (2.0 * P)
        if not c > 0:
            raise ValueError(f"c-hat = {c!r} is not positive; check the SNP filter")

    order = np.argsort(-(m_g * (r_g + 1)))
    chunks = [order[t::n_threads] for t in range(n_threads)]

    def f(B):
        B = np.asarray(B, dtype=np.float64)
        B2 = B.reshape(n, -1)
        k = B2.shape[1]
        B32 = B2.astype(np.float32)

        def part(idx):
            acc = np.zeros((n, k))
            for i in idx:
                cols, ctr, scl, T, lam, _, full = genes[i]
                Z = (G[:, cols].astype(np.float32) - ctr) * scl
                mg = Z.shape[1]
                sq = np.sqrt(lam)
                Qs = (Z @ T) * sq                               # Z V_r
                s_step = max(1, max_cols // mg)
                for a in range(0, Qs.shape[1], s_step):
                    Qa = Qs[:, a:a + s_step]
                    ns = Qa.shape[1]
                    A = (Qa[:, :, None] * Z[:, None, :]).reshape(n, ns * mg)
                    acc += A @ (A.T @ B32)
                if full:
                    D = Z * Z
                else:
                    D = (Qs @ (T * sq).T) * Z                   # Dt = ((Z V_r) V_r') o Z
                acc -= D @ (D.T @ B32)
            return acc

        with ThreadPoolExecutor(n_threads) as ex:
            out = sum(ex.map(part, chunks)) / (2.0 * P * c)
        return out.reshape(B.shape)

    W = spla.LinearOperator((n, n), matvec=f, matmat=f, rmatvec=f, rmatmat=f, dtype=np.float64)
    W.P, W.n_genes, W.rank, W.c = P, len(genes), r_g, c
    return W
def fill_missing_rounded_mean(G, block=4096):
    """In place: each missing call (-9) of G becomes the rounded mean of its SNP.

    int8 cannot hold the exact mean; the rounded mean keeps G integer.  Done
    block by block so the temporary arrays stay small.
    """
    n, m = G.shape
    for a in range(0, m, block):
        g = G[:, a:a + block]
        obs = g >= 0
        mu0 = np.where(obs, g, 0).sum(axis=0, dtype=np.int64) / obs.sum(axis=0)
        i_na, j_na = np.nonzero(~obs)
        g[i_na, j_na] = np.rint(mu0[j_na]).astype(np.int8)
    return G


def genetic_kernels(G, components=('A', 'D'), block=512, n_threads=1):
    """Additive and dominance kernels as LinearOperators sharing ONE pass over G.

    'A': K_A = Z Z' / m,  Z = (G - 2p) / sqrt(2pq)   (HWE standardization, p = mean / 2)
    'D': K_D = D D' / m,  D = GCTA dominance code: g = 0 -> -p/q, 1 -> 1, 2 -> -q/p
    Returns a list of (n, n) LinearOperators, one per component, in the order given.
    All operators share one core pass: applying them to the same B (e.g. inside
    V = sA*KA + sD*KD + se*I) reads G once; later calls with an equal B hit a 1-entry cache.
    Blocks are split over n_threads threads (set BLAS to 1 thread). Monomorphic SNPs get weight 0.

    G: (n, m) int8 with no missing calls (fill_missing_rounded_mean first).
    """
    if not set(components) <= {'A', 'D'} or len(components) == 0:
        raise ValueError(f"components must be a non-empty subset of ('A', 'D'), got {components}")
    n, m = G.shape
    starts = list(range(0, m, block))
    chunks = [starts[t::n_threads] for t in range(n_threads)]

    mu = np.empty(m)

    def stats(blocks):
        for a in blocks:
            mu[a:a + block] = G[:, a:a + block].mean(axis=0, dtype=np.float64)

    with ThreadPoolExecutor(n_threads) as ex:
        list(ex.map(stats, chunks))

    used = (mu > 0) & (mu < 2)
    m_used = int(used.sum())
    p = np.where(used, mu / 2, 0.5)
    q = 1 - p
    w = np.zeros(m)
    w[used] = 1.0 / np.sqrt(2 * p[used] * q[used])
    tab = np.zeros((3, m), dtype=np.float32)
    tab[0, used] = -p[used] / q[used]
    tab[1, used] = 1.0
    tab[2, used] = -q[used] / p[used]

    do_A = 'A' in components
    do_D = 'D' in components

    def core(B2):                                       # (n, k) float64 -> {c: K_c B}
        k = B2.shape[1]
        B32 = B2.astype(np.float32)
        colsum = B2.sum(axis=0)

        def part(blocks):
            acc_A = np.zeros((n, k)) if do_A else None
            corr_A = np.zeros(k) if do_A else None
            acc_D = np.zeros((n, k)) if do_D else None
            d = np.empty((n, block), dtype=np.float32, order='F') if do_D else None
            for a in blocks:
                Gb = G[:, a:a + block]
                nb = Gb.shape[1]
                if do_A:
                    g = Gb.astype(np.float32)
                    wb, mb = w[a:a + block], mu[a:a + block]
                    WC = (g.T @ B32 - mb[:, None] * colsum[None, :]) * (wb * wb)[:, None]
                    acc_A += g @ WC.astype(np.float32)
                    corr_A += mb @ WC
                if do_D:
                    for j in range(nb):
                        d[:, j] = tab[:, a + j][Gb[:, j]]
                    db = d[:, :nb]
                    acc_D += db @ (db.T @ B32)
            return acc_A, corr_A, acc_D

        with ThreadPoolExecutor(n_threads) as ex:
            parts = list(ex.map(part, chunks))
        res = {}
        if do_A:
            res['A'] = (sum(t[0] for t in parts) - sum(t[1] for t in parts)[None, :]) / m_used
        if do_D:
            res['D'] = sum(t[2] for t in parts) / m_used
        return res

    cache = {'B': None, 'res': None}

    def apply(B2):
        hit = (cache['B'] is not None and cache['B'].shape == B2.shape
               and np.array_equal(cache['B'], B2))
        if not hit:
            cache['B'] = B2.copy()
            cache['res'] = core(B2)
        return cache['res']

    def make_op(c):
        def f(B):
            B = np.asarray(B, dtype=np.float64)
            return apply(B.reshape(n, -1))[c].reshape(B.shape).copy()
        return spla.LinearOperator((n, n), matvec=f, matmat=f, rmatvec=f, rmatmat=f, dtype=np.float64)

    return [make_op(c) for c in components]