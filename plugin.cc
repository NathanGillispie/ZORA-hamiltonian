/*
 * @BEGIN LICENSE
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
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/molecule.h"

namespace psi{ namespace zora_core_excitation {

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
SharedWavefunction zora_core_excitation(SharedWavefunction ref_wfn, Options& options)
{
    int print = options.get_int("PRINT");
    printf("Print option = %d\n", print);

    if(!ref_wfn) throw PSIEXCEPTION("SCF has not ran yet!");

    // Need to check DDX options

    //MintsHelper MintsHelper(ref_wfn->basisset(), options, print);

    std::shared_ptr<Molecule> molecule = ref_wfn->molecule();
    std::shared_ptr<BasisSet> aoBasis = ref_wfn->basisset();
    std::shared_ptr<SOBasisSet> soBasis = ref_wfn->sobasisset();

    //molecule->print();

    aoBasis->initialize_singletons();

    const Dimension dimension = soBasis->dimension();

    auto integral = std::make_shared<IntegralFactory>(aoBasis, aoBasis, aoBasis, aoBasis);

    // The matrix factory can create matrices of the correct dimensions...
    auto factory = std::make_shared<MatrixFactory>();
    factory->init_with(dimension, dimension);

    // Form the one-electron integral objects from the integral factory
    std::shared_ptr<OneBodySOInt> sOBI(integral->so_overlap());
    std::shared_ptr<OneBodySOInt> tOBI(integral->so_kinetic());
    std::shared_ptr<OneBodySOInt> vOBI(integral->so_potential());
    // Form the one-electron integral matrices from the matrix factory
    SharedMatrix sMat(factory->create_matrix("Overlap"));
    SharedMatrix tMat(factory->create_matrix("Kinetic"));
    SharedMatrix vMat(factory->create_matrix("Potential"));
    SharedMatrix hMat(factory->create_matrix("One Electron Ints"));

    // Compute the one electron integrals, telling each object where to
    // store the result
    sOBI->compute(sMat);
    tOBI->compute(tMat);
    vOBI->compute(vMat);

    if(print > 5){
        sMat->print();
    }
    if(print > 3){
        tMat->print();
        vMat->print();
    }
    // Form h = T + V by first cloning T and then adding V
    hMat->copy(tMat);
    hMat->add(vMat);
    hMat->print();

    return ref_wfn;
}

}}

