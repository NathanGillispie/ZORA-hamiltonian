#
# @BEGIN LICENSE
#
# zora_core_excitation by Nathan Gillispie, a plugin to:
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

"""
Test plugin initialization file. Read as a python module by psi4.
"""
__version__ = '0.1'
__author__  = 'Nathan Gillispie'

# Load Python modules
from .pymodule import *

# Load C++ plugin
import os
import psi4
plugdir = os.path.split(os.path.abspath(__file__))[0]
sofile = plugdir + '/' + os.path.split(plugdir)[1] + '.so'
psi4.core.plugin_load(sofile)