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
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libfock/points.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/v.h"
#include "psi4/libscf_solver/hf.h"
#include <cmath>

#ifdef _OPENMP
#include "psi4/libpsi4util/process.h"
#include <omp.h>
#endif

namespace psi{ namespace zora_core_excitation {

#include "bigScaryModelBasis.h"
#define C 137.037


void compute_veff(std::shared_ptr<Molecule> mol, std::shared_ptr<DFTGrid> &grid, SharedMatrix veff) {

	int nthreads_ = 1;
#ifdef _OPENMP
	nthreads_ = Process::environment.get_n_threads();
#endif
	int natoms = mol->natom();
	for (int a = 0; a < natoms; a++) {
		int Z = mol->Z(a);
		if (Z > 104) throw PSIEXCEPTION("Z too big. Max value 104");
		auto pos_a = mol->xyz(a);
		
		//Get the list of coefficients and alphas for given atom
		double* coef_a  = &coeffs[c_aIndex[Z-1]];
		double* alpha_a = &alphas[c_aIndex[Z-1]];
		int nc_a = c_aIndex[Z] - c_aIndex[Z-1];

		int index = 0;
#pragma omp parallel for schedule(dynamic) num_threads(nthreads_)
		for (const auto &block : grid->blocks()) {
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
				veff->set(index, p, outer);
			}
			index++;
		}
	}
}

//Scalar Relativistic Kinetic Energy Matrix
void compute_TSR(std::shared_ptr<DFTGrid> &grid, BasisFunctions &props, SharedMatrix veff, SharedMatrix &T_SR) {
	double** T_SRp = T_SR->pointer();
	int max_funcs = props.max_functions();

	double* kernel = new double[props.max_points()];

	int index = 0;
	for (const auto &block : grid->blocks()) {
	  const auto &bf_map = block->functions_local_to_global();
		auto local_nbf = bf_map.size();
		int npoints = block->npoints();

		props.compute_functions(block);
		auto phi_x = props.basis_value("PHI_X");
		auto phi_y = props.basis_value("PHI_Y");
		auto phi_z = props.basis_value("PHI_Z");

		double* w = block->w();

		//preprocess kernel c²/(2c²-veff) * weight
		for (int p = 0; p < npoints; p++) {
			kernel[p] = C *C /(2.*C *C - veff->get(index,p)) * w[p];
		}

		for (int l_mu = 0; l_mu < local_nbf; l_mu++) {
			int mu = bf_map[l_mu];
			for (int l_nu = l_mu; l_nu < local_nbf; l_nu++) {
				int nu = bf_map[l_nu];
				for (int p = 0; p < npoints; p++) {
					T_SRp[mu][nu] += ( phi_x->get(p,l_mu)*phi_x->get(p,l_nu) + phi_y->get(p,l_mu)*phi_y->get(p,l_nu) + phi_z->get(p,l_mu)*phi_z->get(p,l_nu) ) * kernel[p];
				}
			}
		}
		index++;
	}

	T_SR->copy_upper_to_lower();
	delete kernel;
}

void compute_SO(std::shared_ptr<DFTGrid> &grid, BasisFunctions &props, SharedMatrix veff, SharedMatrix H_SOx, SharedMatrix H_SOy, SharedMatrix H_SOz) {
	int max_funcs = props.max_functions();
	double** H_SOxp = H_SOx->pointer();
	double** H_SOyp = H_SOy->pointer();
	double** H_SOzp = H_SOz->pointer();

	double* kernel = new double[props.max_points()];

	int index = 0;
	for (const auto &block : grid->blocks()) {
		int npoints = block->npoints();
		int local_nbf = block->local_nbf();

		props.compute_functions(block);
		auto phi_x = props.basis_value("PHI_X");
		auto phi_y = props.basis_value("PHI_Y");
		auto phi_z = props.basis_value("PHI_Z");

		auto w = block->w();
		//preprocess kernel veff/(4c²-2veff) * weight
		for (int p = 0; p < npoints; p++) {
			kernel[p] = veff->get(index, p)/(4.*C *C  - 2.*veff->get(index, p)) * w[p];
		}
		
	  const auto &bf_map = block->functions_local_to_global();
		//compute upper triangular half. H_SO are antisymmetric.
		for (int mu_l = 0; mu_l < local_nbf; mu_l++) {
			int mu = bf_map[mu_l];
			for (int nu_l = mu+1; nu_l < local_nbf; nu_l++) {
				int nu = bf_map[nu_l];
				for (int p = 0; p < npoints; p++) {
					H_SOxp[mu][nu] += (phi_y->get(p,mu_l)*phi_z->get(p,nu_l) - phi_z->get(p,mu_l)*phi_y->get(p,nu_l)) * kernel[p];
					H_SOyp[mu][nu] += (phi_z->get(p,mu_l)*phi_x->get(p,nu_l) - phi_x->get(p,mu_l)*phi_z->get(p,nu_l)) * kernel[p];
					H_SOzp[mu][nu] += (phi_x->get(p,mu_l)*phi_y->get(p,nu_l) - phi_y->get(p,mu_l)*phi_x->get(p,nu_l)) * kernel[p];
				}
			}
		}
		index++;
	}
	
	//keep it simple silly
	for (int mu = 0; mu < max_funcs; mu++) {
		for (int nu = mu - 1; nu >= 0; nu--) {
			H_SOxp[mu][nu] = -H_SOxp[nu][mu];
			H_SOyp[mu][nu] = -H_SOyp[nu][mu];
			H_SOzp[mu][nu] = -H_SOzp[nu][mu];
		}
	}
	
	delete kernel;
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
	auto Vpot = ref_wfn->V_potential();
	if (!Vpot) throw PSIEXCEPTION("Must run DFT method");

	std::shared_ptr<DFTGrid> grid = Vpot->grid();
	//would ideally like to replace this grid object

	auto primary = ref_wfn->basisset();
	BasisFunctions bf_computer(primary, grid->max_points(), grid->max_functions());
	bf_computer.set_deriv(1);
	int max_funcs = grid->max_functions();
	std::cout << "max functions: " << max_funcs << std::endl;

	timer_on("Compute Veff");
	int nblocks = grid->blocks().size();
	auto veff = std::make_shared<Matrix>("Effective potential", nblocks, bf_computer.max_points());
	compute_veff(primary->molecule(), grid, veff);
	timer_off("Compute Veff");

	printf("Computed effective potential.\n");

	timer_on("Scalar Relativistic Kinetic");
	auto T_SR = std::make_shared<Matrix>(max_funcs, max_funcs);
	compute_TSR(grid, bf_computer, veff, T_SR);
	timer_off("Scalar Relativistic Kinetic");

	printf("Computed TSR.\n");

	timer_on("Spin-orbit terms");
	auto H_SOx = std::make_shared<Matrix>("H_SOx", max_funcs, max_funcs);
	auto H_SOy = std::make_shared<Matrix>("H_SOy", max_funcs, max_funcs);
	auto H_SOz = std::make_shared<Matrix>("H_SOz", max_funcs, max_funcs);
	compute_SO(grid, bf_computer, veff, H_SOx, H_SOy, H_SOz);
	timer_off("Spin-orbit terms");

	printf("Computed spin-orbit terms.\n");

	T_SR->save("mats/T_SR.mat", false, false);
	H_SOx->save("mats/H_SOx.mat", false, false);
	H_SOy->save("mats/H_SOy.mat", false, false);
	H_SOz->save("mats/H_SOz.mat", false, false);

	return ref_wfn;
}

}}

