from Function_input import *
from Function_estimation import mc_reml, reml_se, h2_se
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--trait', type=str, required=True)
# Chromosomes, e.g. --chroms 22  or  --chroms $(seq 1 22)
parser.add_argument('--chroms', type=int, nargs='+', default=[22])
# Genetic components to estimate (E is always included), e.g. --components A D GxG
parser.add_argument('--components', type=str, nargs='+', default=['A'],
                    choices=['A', 'D', 'GxG'])
# Paths
parser.add_argument('--geno_dir', type=str,
                    default='/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_100K_HM3')
parser.add_argument('--phen_dir', type=str,
                    default='/dbgap/martinzh_ukb16549/ziyanzha/Phenotype/phen')
parser.add_argument('--cov_dir', type=str,
                    default='/dbgap/martinzh_ukb16549/martin/UKBB_trait_from_Steven')
parser.add_argument('--gene_dir', type=str,                 # long-format gene tables, one per chromosome (GxG)
                    default='/dbgap/martinzh_ukb16549/ziyanzha/Genotype/UKBimp_337K_HM3/gxg_gene_table')
parser.add_argument('--out_dir', type=str,
                    default='/dbgap/martinzh_ukb16549/ziyanzha/results/MCREML_real')
# GxG kernel
parser.add_argument('--gxg_rank', type=int, default=9999)     # rank r of the GxG kernel; 9999 = exact
# MC-REML
parser.add_argument('--nmc', type=int, default=100)          # Hutchinson probes (fine phase)
parser.add_argument('--nmc_coarse', type=int, default=15)    # probes in the coarse phase
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--verbose', action='store_true')        # print the REML iterations
parser.add_argument('--threads', type=int, default=1)        # threads for the kernel products (set BLAS to 1 thread)
args = parser.parse_args()

trait = args.trait
components = [c for c in ('A', 'D', 'GxG') if c in args.components]   # fixed order: A, D, GxG
if 'GxG' in components and args.gene_dir is None:
    parser.error('--gene_dir is required for GxG')

pgen_prefixes = ['%s/ukb_imp_chr%d_v3' % (args.geno_dir, c) for c in args.chroms]
phen_file = '%s/%s.phen' % (args.phen_dir, trait)
if trait.startswith('biochemistry_'):
    cov_file = args.cov_dir + '/UKB_biochemistry.051619.cov'
else:
    cov_file = args.cov_dir + '/ukb4777.processed_and_post.plinkPCs.tab.gz'

# ---------------- samples, phenotype, covariates ----------------
iids = read_psam_iids(pgen_prefixes[0])
y_all = read_phenotype(phen_file, trait, iids)
df_reg = build_covariate_table(cov_file, trait)
rows, y, X, x_cols = select_samples(iids, y_all, df_reg)

# ---------------- genotype: int8, standardized on the fly ----------------
# MAF > 1% is applied upstream (plink2), so no MAF filter here
G, snp_chrom, snp_idx = read_genotype_chroms(pgen_prefixes, rows, iids)
fill_missing_rounded_mean(G)

# ---------------- kernels ----------------
# A and D in one pass over G
comps_AD = [c for c in components if c != 'GxG']
K_AD = genetic_kernels(G, comps_AD, n_threads=args.threads) if comps_AD else None

# Within-gene AxA; gene sets index the columns of G (unfiltered snp_chrom / snp_idx)
K_gxg = None
if 'GxG' in components:
    gene_files = ['%s/gene_table_chr%d.tsv' % (args.gene_dir, c) for c in args.chroms]
    gene_ids, gene_sets = read_gene_sets(gene_files, pgen_prefixes, snp_chrom, snp_idx)
    K_gxg = gxg_kernel(G, gene_sets, r=args.gxg_rank, n_threads=args.threads)


def K_all(B):
    """[K_A B, K_D B, K_GxG B] for the requested components, in that order."""
    out = K_AD(B) if K_AD is not None else []
    if K_gxg is not None:
        out.append(K_gxg(B))
    return out


# ---------------- MC-REML: V = sum_c s2c K_c + s2e I, covariates X as fixed effects ----------------
s, AI, info = mc_reml(K_all, y, X, Nmc=args.nmc, Nmc_coarse=args.nmc_coarse,
                      seed=args.seed, verbose=args.verbose)
se = reml_se(AI, args.nmc, info['at_bound'])
h2, se_h2 = h2_se(s, AI, args.nmc, info['at_bound'])

# ---------------- output ----------------
chrom_label = ('chr%d' % args.chroms[0] if len(args.chroms) == 1
               else 'chr%d-%d' % (args.chroms[0], args.chroms[-1]))
os.makedirs(args.out_dir, exist_ok=True)
comp_label = ''.join(components)                           # e.g. A, AD, ADGxG
out_path = '%s/%s_%s_%s.tsv' % (args.out_dir, trait, chrom_label, comp_label)
row = {'trait': trait, 'chroms': chrom_label, 'components': comp_label,
       'n': len(y), 'm': G.shape[1], 'n_covar': X.shape[1]}
for i, c in enumerate(components):                         # s order: components, then e
    row['s2' + c.lower()] = s[i]
    row['se_s2' + c.lower()] = se[i]
row['s2e'] = s[-1]
row['se_s2e'] = se[-1]
for i, c in enumerate(components):
    row['h2' + c.lower()] = h2[i]
    row['se_h2' + c.lower()] = se_h2[i]
if K_gxg is not None:
    row.update({'gxg_n_genes': K_gxg.n_genes, 'gxg_P': K_gxg.P, 'gxg_c': K_gxg.c,
                'gxg_r': args.gxg_rank, 'gxg_max_rank': int(K_gxg.rank.max())})
row.update({'nmc': args.nmc, 'converged': info['converged'],
            'n_iters': info['n_iters'], 'n_reject': info['n_reject'],
            'at_bound': ','.join(str(int(b)) for b in info['at_bound'])})
pd.DataFrame([row]).to_csv(out_path, sep='\t', index=False)