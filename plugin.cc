/* @BEGIN LICENSE
 *
 * zora_core_excitation by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/matrix.h"

#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libfock/cubature.h"

namespace psi{ namespace zora_core_excitation {

void compute_veff(std::shared_ptr<Molecule> mol, std::shared_ptr<DFTGrid> grid, SharedMatrix ret) {
	int natoms = mol->natom();
	//helper->potential_integral(std::vector<SharedVector>& grid data)
	//molecule->Z(int atom)
	//molecule->xyz(int atom) -> Vector3
}

extern "C" PSI_API
int read_options(std::string name, Options& options)
{
    //Called when input file has a set mymodule key value
    //Comments must be added in /*- -*/ form for documentation
    if (name == "ZORA_CORE_EXCITATION"|| options.read_globals()) {
        /*- Used to specify how much is printed to the output -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" PSI_API
SharedWavefunction zora_core_excitation(std::shared_ptr<scf::HF> ref_wfn, Options& options) {
	if(!ref_wfn) throw PSIEXCEPTION("SCF has not ran yet!");

	int print = options.get_int("PRINT");
	printf("Print option = %d\n", print);

	auto Vpot = ref_wfn->V_potential();
	if (!Vpot) throw PSIEXCEPTION("Must run DFT method");
	Vpot->initialize();
	auto grid = Vpot->grid();

	std::shared_ptr<PointFunctions> props = Vpot->properties()[0];
	props->set_pointers(ref_wfn->Da());

	auto veff = std::make_shared<Matrix>(Vpot->nblocks(), grid->max_points(), 3);
	compute_veff(ref_wfn->molecule(), grid, veff);

	//Compute integrals, requires evaluating phi_mu(r)
	//Can be done with the numinthelper line but I'm starting to think
	//I would rather just do the integration myself. There's no reason to
	//use std::vector for this. I know exactly what size and datatype I need

	Vpot->finalize();
	return ref_wfn;
}

}}

