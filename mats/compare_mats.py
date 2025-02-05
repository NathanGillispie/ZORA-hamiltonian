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

    bohr_keys = {'T_SR', 'H_SOx', 'H_SOy', 'H_SOz'}
    common = list(set(psi4m.keys()).intersection(bohr_keys))

    if (bohrm[0].shape != psi4m[common[0]].shape):
        exit("Shape mismatch")

    for i, name in enumerate(common):
        print ("\n  Diffs of", name)
        diff_mat(psi4m[name], bohrm[i])
    np.set_printoptions(linewidth = 100, floatmode="fixed", precision=6,suppress=True)

    a = psi4m['T_SR']
    b = bohrm[0]

    print("psi4m\n",np.array_str(a))
    #print("bohrm\n",np.array_str(b))
