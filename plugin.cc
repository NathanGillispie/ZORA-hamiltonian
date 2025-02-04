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
#include "psi4/libmints/vector.h"

#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libscf_solver/hf.h"
#include "psi4/libfock/cubature.h"
#include <cmath>
#include <vector>

namespace psi{ namespace zora_core_excitation {

#include "bigScaryModelBasis.h"

void compute_veff(std::shared_ptr<Molecule> mol, std::shared_ptr<VBase> Vpot, SharedMatrix ret) {
	double** veff = ret->pointer();
	int nblocks = ret->nrow();
	int max_pts = ret->ncol();
	
	int natoms = mol->natom();
	for (int a = 0; a < natoms; a++) {
		int Z = mol->Z(a);
		if (Z > 104) throw PSIEXCEPTION("Choose a chemically relevant system (Z too big)");
		auto pos_a = mol->xyz(a);

		double* coef_a  = &coeffs[c_aIndex[Z-1]];
		double* alpha_a = &alphas[c_aIndex[Z-1]];
		int nc_a = c_aIndex[Z] - c_aIndex[Z-1];

	//pragma parallelize stuff
		for (int b = 0; b < nblocks; b++) {
			auto block = Vpot->get_block(b);
			int npoints = block->npoints();
			double* x = block->x();
			double* y = block->y();
			double* z = block->z();

			//einsums("i,ip->p", 𝕔[A], erf(α⊗ r))/r
			for (int p = 0; p < npoints; p++) {
				double dist = hypot(pos_a[0]-x[p], pos_a[1]-y[p], pos_a[2]-z[p]);
				double outer = 0;
				for (int i = 0; i < nc_a; i++) {
					//check to make sure this is the correct erf
					outer += std::erf(dist * alpha_a[i]) * coef_a[i];
				}
				outer /= dist;
				outer -= Z/dist;
				veff[b][p] = outer;
			}
	//pragma parallelize
		}
	}
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

	std::shared_ptr<PointFunctions> props = Vpot->properties()[0];
	//may not be necessary to set Da pointer.
	props->set_pointers(ref_wfn->Da());
	props->set_ansatz(0); //necessary to avoid segfault when computing points
	props->set_deriv(1);
	
	timer_on("Compute Veff");
	auto veff = std::make_shared<Matrix>(Vpot->nblocks(), props->max_points());
	compute_veff(ref_wfn->molecule(), Vpot, veff);
	timer_off("Compute Veff");
	
	int max_funcs = props->max_functions();
	auto T_SR = std::make_shared<Matrix>(max_funcs, max_funcs);

	timer_on("Compute Scalar Relativistic T");
	double** T_SRp = T_SR->pointer();
	double**  veffp = veff->pointer();
	double kernel[props->max_points()];

#define C 137.037
	for (int b = 0; b < Vpot->nblocks(); b++) {
		auto block = Vpot->get_block(b);
		int npoints = block->npoints();

		props->compute_points(block);
		double** phi_x = props->basis_value("PHI_X")->pointer();
		double** phi_y = props->basis_value("PHI_Y")->pointer();
		double** phi_z = props->basis_value("PHI_Z")->pointer();

		auto w = block->w();

		//preprocess kernel c^2/(2c^2-veff) * weight
		for (int p = 0; p < npoints; p++) {
			kernel[p] = C*C/(2.*C*C - veffp[b][p]) * w[p];
		}

		for (int mu = 0; mu < max_funcs; mu++) {
			for (int nu = 0; nu < max_funcs; nu++) {
				for (int p = 0; p < npoints; p++) {
					T_SRp[mu][nu]+= (phi_x[p][mu] * phi_x[p][nu] + phi_y[p][mu] * phi_y[p][nu] + phi_z[p][mu] * phi_z[p][nu]) * kernel[p];
				}
			}
		}
	}
	timer_off("Compute Scalar Relativistic T");

	Vpot->finalize();
	
	T_SR->print_out();

	return ref_wfn;
}

}}

