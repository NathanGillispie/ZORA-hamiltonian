#!/usr/bin/env python

"""
Psi4 + numpy test to get stuff working
"""

import psi4
import numpy as np
import os
from scipy import special

# Some valid spherical points options:
# 50, 86, 110, 194, 230, 350

OPTIONS = {
        'basis': 'sto-3g',
        'scf_type': 'pk',
        'reference': 'rks',
        'e_convergence': 1e-6,
        'dft_spherical_points': 194,
        'dft_radial_points': 200,
        'd_convergence': 1e-6
}

FAILING_MOL = """
3 2
V  0 0 0 
F  1.86 0 0
F  -.93 1.61 0
F   .93 -1.61 0
symmetry c1
units angstrom
"""

WATER="""
0 1
  O      0.00000000     0.00000000     0.11726921
  H      0.75698224     0.00000000    -0.46907685
  H     -0.75698224     0.00000000    -0.46907685
units angstrom
symmetry c1
"""

def getDeriv(Vpot, grid):
    max_funcs = Vpot.properties()[0].max_functions()
    deriv = np.empty((0,4,max_funcs))
    bf_computer = Vpot.properties()[0]
    for block in grid.blocks():
        bf_computer.compute_functions(block)

        npts = block.npoints()
        phi   = bf_computer.basis_values()['PHI'  ].to_array()[:npts]
        phi_x = bf_computer.basis_values()['PHI_X'].to_array()[:npts]
        phi_y = bf_computer.basis_values()['PHI_Y'].to_array()[:npts]
        phi_z = bf_computer.basis_values()['PHI_Z'].to_array()[:npts]

        ex = np.array((phi,phi_x,phi_y,phi_z)).transpose((1,0,2))
        deriv = np.concat((deriv,ex))
    deriv = deriv.transpose((1,0,2))
    return (deriv[0], deriv[1], deriv[2], deriv[3])

def getWeights(Vpot, grid):
    max_funcs = Vpot.properties()[0].max_functions()
    weights = np.empty((0))
    bf_computer = Vpot.properties()[0]
    for block in grid.blocks():
        npts = block.npoints()
        w = block.w().to_array()[:npts]

        weights = np.concat((weights,w))
    return weights

def init(mol_str):
    psi4.set_memory(int(2e9))
    psi4.core.set_output_file('output.dat', False)

    psi4.set_options(OPTIONS)

    molecule = psi4.geometry(mol_str)

    psi4.core.Wavefunction.build(molecule, psi4.core.get_global_option('basis'))
    e, wfn = psi4.energy('pbe0', return_wfn=True)
    print("  pbe0 energy:",e)
    return wfn

def read_basis(atoms):
    '''
    Reads the two-component model basis file used to generate Veff.
    Returns coefficients and sqrt alphas
    '''
    basis_file = []
    with open(os.path.abspath(os.path.dirname(__file__)+"modbas.2c"),'r') as f:
      for line in f.readlines():
        if (len(line) > 1):
          basis_file.append(line)

    c_a = []
    for atom in map(lambda a: a[0].lower(), atoms):
      # get position of given atom in basis_file
      position = [line for line, a in enumerate(basis_file) if a.split()[0] == atom][0]
      self.nbasis = int(basis_file[position][10:15])
      # Assume the basis_file has the same delimiters
      array = np.loadtxt(basis_file[position+2:position+2+self.nbasis]).transpose((1,0))
      c_a.append((np.array(array[1]),np.sqrt(np.array(array[0]))))
    return c_a

def sameMat(a, b):
    diff = a-b
    e = np.sum(diff.flat)
    #arbitrary cutoff
    return bool(e < np.float64(1e-7))

def compute_veff(mol, grid):
    nblocks = 0
    npoints = grid.maxpoints()
    veff = np.zeros((nblocks, points))

    natoms = mol.natom()
    atom_crds = mol.to_arrays() # mol.xyz(atom)
    atom_strings = atom_crds[2].tolist()
    atom_charges = atom_crds[3].tolist()

    xyzw = Vpot.get_np_xyzw()
    w = xyzw[3]
    npts = w.shape[0]
    xyz = np.pose(np.array(xyzw[:-1]))
    RPA = np.sqrt(np.sum(xyz**2, axis=1))

    return veff

wfn = init(WATER)

mol = wfn.molecule()
basis = wfn.basisset()

Vpot = wfn.V_potential()
grid = Vpot.grid().build(mol, basis)

# bf_computer = Vpot.properties()[0]
bf_computer = psi4.core.BasisFunctions(basis, grid.max_points(), grid.max_functions())

max_funcs = basis.nbf()

integral = np.zeros((max_funcs,max_funcs))

"""Compute the total density"""
for block in grid.blocks():
    bf_map = block.functions_local_to_global()
    local_nbf = len(bf_map)
    local_npts = block.npoints()

    w = block.w().to_array()
    bf_computer.compute_functions(block)
    phi = bf_computer.basis_values()['PHI'].to_array()[:local_npts]

    for l_mu in range(local_nbf):
        mu = bf_map[l_mu]
        for l_nu in range(local_nbf):
            nu = bf_map[l_nu]
            for p in range(local_npts):
                integral[mu][nu] += w[p] * phi[p][l_mu] * phi[p][l_nu]


print_opts = {
        "max_line_width": 200,
        "floatmode": "fixed",
        "precision": 4,
        "suppress_small" : True
}

print(np.array2string(integral, **print_opts))
print()
print(np.array2string(np.diag(integral), **print_opts))

