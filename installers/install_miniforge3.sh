#!/bin/bash

# WildSEQ: a qPCR Sequencing Bioinformatics Pipeline
# Author: S. Sturrock
#
#
# This file is part of WildSEQ.
# Copyright (C) 2025 Shane Sturrock
#
# WildSEQ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


# Script to auto install Miniforge3 for users who don't already have it

ROOT=${HOME}
CONDA_VER=25.3.0-3
CONDA_BIN=${ROOT}/miniforge3/bin

# Delete old miniforge3 if it exists
if [ -d ${ROOT}/miniforge3 ]; then
	rm -rf ${ROOT}/miniforge3
fi

# Download and install miniforge
wget -N https://github.com/conda-forge/miniforge/releases/download/${CONDA_VER}/Miniforge3-${CONDA_VER}-Linux-x86_64.sh
bash ./Miniforge3-${CONDA_VER}-Linux-x86_64.sh -b -p ${ROOT}/miniforge3/

# Install the conda setup in user's .bashrc
${CONDA_BIN}/conda init

# Clean up
rm Miniforge3-${CONDA_VER}-Linux-x86_64.sh
