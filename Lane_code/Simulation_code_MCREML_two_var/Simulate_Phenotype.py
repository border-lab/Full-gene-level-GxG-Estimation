from Function_MCREML import *
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--m', type=int, required=True)
parser.add_argument('--n', type=int, required=True)
parser.add_argument('--s2a', type=float, required=True)
parser.add_argument('--s2d', type=float, required=True)
parser.add_argument('--s2e', type=float, required=True)
parser.add_argument('--rep', type=int, required=True)
parser.add_argument('--mode', type=str, required=True)


args = parser.parse_args()

m = args.m
n = args.n
s2a = args.s2a
s2d = args.s2d
s2e = args.s2e
mode = args.mode
rep = args.rep


# Read genotype
SNP = pd.read_csv(f"/home/ziyanzha/MOM_within_gene/stored_genotype/{mode}SNP_n{n}_m{m}.csv", header=None)
SNP = SNP.to_numpy()

# Load the precomputed Cholesky factors (La La' = s2a K_a, Ld Ld' = s2d K_d).
save_dir = "/home/ziyanzha/MOM_within_gene/MCREML_two_var/Cholesky"
La = np.load(f"{save_dir}/La_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2e{s2e}.npy")
Ld = np.load(f"{save_dir}/Ld_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2e{s2e}.npy")

# Additive + dominance phenotype: y = g_a + g_d + e with exact-variance rescaling.
Za, Zd, y = simulate_remove_sampling_err(SNP, La, Ld, s2a=s2a, s2d=s2d, s2e=s2e)

output_dir = f"/home/ziyanzha/MOM_within_gene/MCREML_two_var/Phenotype/y_{mode}_n{n}_m{m}_s2a{s2a}_s2d{s2d}_s2e{s2e}"

os.makedirs(output_dir, exist_ok=True)
y_path = f"{output_dir}/rep{rep}.csv"
pd.DataFrame(y).to_csv(y_path, index=False, header=False)


# Save Za, Zd only once (both are deterministic for a given genotype matrix).
z_dir = "/home/ziyanzha/MOM_within_gene/stored_genotype"
za_path = f"{z_dir}/Za_{mode}_n{n}_m{m}.csv"
zd_path = f"{z_dir}/Zd_{mode}_n{n}_m{m}.csv"
os.makedirs(z_dir, exist_ok=True)
if not os.path.exists(za_path):
    pd.DataFrame(Za).to_csv(za_path, index=False, header=False)
    print(f"Za saved to: {za_path}")
if not os.path.exists(zd_path):
    pd.DataFrame(Zd).to_csv(zd_path, index=False, header=False)
    print(f"Zd saved to: {zd_path}")
