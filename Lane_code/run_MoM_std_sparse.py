from function import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2gxg', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--nmc', type=int, default=40)
parser.add_argument('--density', type=float, default=0.1)  # Fixed: removed space
parser.add_argument('--rep', type=int, required=True)
args = parser.parse_args()

m = args.m
n = args.n
s2gxg = args.s2gxg
s2e = args.s2e
nmc = args.nmc
rep = args.rep
density = args.density

Random1KSNP= pd.read_csv("data/chr1_Random1KSNP.raw", sep=r'\s+')
real_data = Random1KSNP.iloc[:, 6:].to_numpy()
real_data = real_data[:n, :m]

Lgxg, Le = simulate_Cholesky_from_std_sparse(real_data, s2gxg, s2e, stability=1e-8, density=density)  
Z, y = simulate_from_raw_only_W_remove_sampling_err(real_data, Lgxg, Le, s2gxg, s2e)
gxg, e, _ = MoM_only_M_std(Z, y, nmc=nmc)

filename = f"Random_n{n}_m{m}_Density{density}.txt"  # Include density in filename
with open(filename, 'a') as f:
    f.write(f"({gxg.item()},{e.item()})\n")