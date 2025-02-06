#!/home/machdev/miniconda3/envs/p4dev/bin/python
import os
import numpy as np
import pickle as pkl

def get_psi4_mats():
    mats_files = [f for f in os.listdir() if f.endswith('.mat')]
    mats = {}

    for file in mats_files:
        nel = -1
        with open(file, 'r') as f:
            f.readline()
            f.readline()
            nel = int(f.readline())
        if (nel <= 0):
            print("Warning:", file, "contains no elements")
            continue
        mat = np.transpose(np.genfromtxt(file, skip_header=3))

        rows = [int(i) for i in mat[0]]
        cols = [int(i) for i in mat[1]]
        dim = max(max(rows), max(cols))

        a = np.zeros((dim+1, dim+1))
        for g in range(mat.shape[1]):
            a[int(mat[0][g])][int(mat[1][g])] = mat[2][g]
        mats[file[:-4]] = a
    return mats

def get_pkl(file):
    with open(file, 'rb') as f:
        return pkl.load(f)
    pass

def save_pkl(file):
    with open(file, 'wb') as f:
        pkl.dump(file, f)

def diff_mat(a, b):
    diff = a - b
    print("Trace of diff:", np.trace(diff))
    print("Abs error of all diffs:", np.sum(np.abs(diff)))
    return diff

if __name__ == '__main__':
    psi4m = get_psi4_mats()
    bohrm = get_pkl('bohrmats.pkl')

    common = list(set(psi4m.keys()).intersection(set(bohrm.keys())))

    if (bohrm[common[0]].shape != psi4m[common[0]].shape):
        print("Shape mismatch:", common[0], "size")
        print("Bohr ->", bohrm[common[0]].shape)
        print("Psi4 ->", psi4m[common[0]].shape)
        exit(1)

    for i, name in enumerate(common):
        print ("\n  Diffs of", name)
        diff_mat(psi4m[name], bohrm[name])
    np.set_printoptions(linewidth = 200, floatmode="fixed", precision=4,suppress=True)

    a = psi4m['T_SR']
    b = bohrm['T_SR']

    print("psi4m\n",np.array_str(a))
    print("bohrm\n",np.array_str(b))
