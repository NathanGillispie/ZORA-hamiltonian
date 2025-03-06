#!/usr/bin/env python

import bohr
import relativistic
from pyscf import dft
import os
import numpy as np

"""RUNS BOHR CALCULATION"""

# find basis sets at
# /home/machdev/miniconda3/envs/bohr/lib/python3.13/site-packages/pyscf/gto/basis


INPUT = """
start molecule
  O      0.00000000     0.00000000     0.11726921
  H      0.75698224     0.00000000    -0.46907685
  H     -0.75698224     0.00000000    -0.46907685
end molecule

charge 0
spin 0
basis def2-tzvp
diis true
units angstrom
relativistic zora
grid_level 8
e_convergence 1e-6
d_convergence 1e-6

xc pbe0
method rks
""".split('\n')


def diff_mat(a, b, ovlp):
    diff = np.abs(a - b)
    print("  Tr : %0.7f"%np.trace(diff))
    print("  Abs: %0.7f"%np.sum(diff))
    print("  Rel: %.2f"%((np.sum(diff)-ovlp)/ovlp))
    return diff

class cringe_rks():
    def __init__(self,mol,options):
        self.ints_factory = mol
        self.options = options
        self.nbf = mol.nao

    def compute(self):
        self.options.xctype  = dft.libxc.xc_type(self.options.xc)
        self.options.xcalpha = dft.libxc.hybrid_coeff(self.options.xc)
        self.jk    = dft.RKS(self.ints_factory, self.options.xc) #jk object from pyscf
        zora = relativistic.ZORA(self)
        zora.get_zora_correction()
        H = {
             "T_SR": zora.T[0],
             "H_SOx": zora.T[1],
             "H_SOy": zora.T[2],
             "H_SOz": zora.T[3],
             "overlap": zora._overlap,
             "T": zora._T,
             "veff": zora._veff
        }
        mol = self.ints_factory

        labels = mol.ao_labels()
        angmom = [int(mol.bas_angular(i)) for i in range(mol.nbas)]

        ao_ang = []
        for l in angmom:
            # for spherical
            nsub = 2*l + 1
            # for cartesian
            # nsub = (l+1)*(l+2)//2
            for _ in range(nsub):
                ao_ang.append(l)
        return H, labels, ao_ang

import pickle
try:
    with open("bohrmats.pkl", "rb") as f:
        bohrm, labels, ao_ang = pickle.load(f)
    print("Successfully loaded pickle. Skipping bohr.")
except FileNotFoundError as err:
    print("Running bohr...")
    opts = bohr.build_options(INPUT)
    mol = bohr.build_mol(opts)

    _wfn = cringe_rks(mol, opts)
    bohrm, labels, ao_ang = _wfn.compute()

    with open("bohrmats.pkl", "wb") as f:
        pickle.dump((bohrm, labels, ao_ang), f)

with open("p4mats.pkl", 'rb') as f:
    psi4m = pickle.load(f)
    print("Successfully loaded psi4 pickle")


np.set_printoptions(precision=4, linewidth=300, suppress=True)

### PRINT T_SR MATS ###
a = psi4m['T_SR']
b = bohrm['T_SR']
print()
if (len(a) > 10):
    print("psi4m (first 10 elements)\n", a[:10,:10])
    print("\nbohrm (first 10 elements)\n", b[:10,:10])
else:
    print("psi4m\n", a)
    print("\nbohrm\n", b)

### GET OVERLAPS ### 
print("\n")
povlp=psi4m['overlap']
bovlp=bohrm['overlap']
nbf = len(ao_ang)

assert nbf == povlp.shape[0]
assert nbf == bovlp.shape[1]

pstart=[]
dstart=[]
fstart=[]
i = 0
while i < nbf:
    L = ao_ang[i]
    if L >= 4:
        exit("G+ orbital permutations unknown, choose a smaller basis.")
    if L == 1:
        pstart.append(i)
        i += 3
        continue
    if L == 2:
        dstart.append(i)
        i += 5
        continue
    if L == 3:
        fstart.append(i)
        i += 7
        continue
    i += 1

# Permutation matricies for aligning psi4 matricies
# For B = bohr matrix,
#     P = psi4 matrix,
# and T = permutation matrix below,
# B = T @ P @ np.transpose(T)
RP = np.array(
        [[0, 1, 0],
        [0, 0, 1],
        [1, 0, 0]], int
)

RD = np.array(
       [[0, 0, 0, 0, 1],
        [0, 0, 1, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 0, 1, 0]], int
)

RF = np.array(
       [[0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0]], int
)

# Total permutation matrix
R = np.zeros((nbf,nbf))
i = 0
while i < nbf:
    L = ao_ang[i]
    if L == 0:
        R[i,i]=1
        i += 1
        continue
    if L == 1:
        for j in range(3):
            for k in range(3):
                R[i+j,i+k] = RP[j,k]
        i += 3
        continue
    if L == 2:
        for j in range(5):
            for k in range(5):
                R[i+j,i+k] = RD[j,k]
        i += 5
        continue
    if L == 3:
        for j in range(7):
            for k in range(7):
                R[i+j,i+k] = RF[j,k]
        i += 7
        continue

assert np.allclose(np.linalg.inv(R), R.T)

ovlperr = np.abs((R@povlp@R.T)-bovlp).sum()
print("Overlap abs error\n  Without perm: %.7f"%np.abs(povlp-bovlp).sum())
print("  With perm   : %.7f"%ovlperr)

for i, name in enumerate(["T_SR", "H_SOx", "H_SOy", "H_SOz"]):
    print(f"\n{name}")
    diff_mat(R @ psi4m[name] @ R.T, bohrm[name], ovlperr)

tdiff = bohrm['T_SR'] - R @ psi4m['T_SR'] @ R.T
for i in range(nbf):
    print(np.array2string(np.diag(tdiff)[i],precision=7),"\t", labels[i])

veffb = bohrm['veff']
veffb = np.asarray([i for i in veffb if i < -1e-9])
veffbmag = np.log10(-1 * veffb)

veff = psi4m['veff'].flatten()
veff = np.asarray([i for i in veff if i < -1e-9])
veffmag = np.log10(-1 * veff)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,2,tight_layout=True)
fig.title = "Histogram of veff"
ax[0].hist(veffmag, bins=20, range=(-5,5))
ax[0].set_xlabel("log10(-veff)")
ax[0].set_title("Psi4")

ax[1].hist(veffbmag, bins=20, range=(-5,5))
ax[1].set_xlabel("log10(-veff)")
ax[1].set_title("Bohr")
plt.show()



