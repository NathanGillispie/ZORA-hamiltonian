#
# @BEGIN LICENSE
#
# zora_core_excitation by Psi4 Developer, a plugin to:
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2023 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

import psi4
import psi4.driver.p4util as p4util
from psi4.driver.procrouting import proc_util

def run_zora_core_excitation(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    zora_core_excitation can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('zora_core_excitation')"""
 
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Compute a SCF reference, a wavefunction is returned.
    # Holds molecule, orbitals, Fock matrices, and more
    ref_wfn = kwargs.get('ref_wfn', None)
    if ref_wfn is None:
        ref_wfn = psi4.driver.scf_helper(name, **kwargs)

    # Ensure IWL files have been written when not using DF/CD
    proc_util.check_iwl_file_from_scf_type(psi4.core.get_option('SCF', 'SCF_TYPE'), ref_wfn)

    psi4.core.timer_on("ZORA")
    zora_wfn = psi4.core.plugin('zora_core_excitation.so', ref_wfn) # Calling plugin: Setting the reference wavefunction in this way is ONLY for plugins
    psi4.core.timer_off("ZORA")

    return zora_wfn

psi4.driver.procedures['energy']['zora'] = run_zora_core_excitation

def form_H(wfn):
    print("potential: ", wfn.V_potential())

    psi4.core.timer_on("ZORA")
    zora_wfn = psi4.core.plugin('zora_core_excitation.so', wfn) # Calling plugin: Setting the reference wavefunction in this way is ONLY for plugins
    psi4.core.timer_off("ZORA")
    return zora_wfn

