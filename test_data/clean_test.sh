#!/usr/bin/env bash

# LAMP Sequencing Bioinformatics Pipeline
# Author: S. Sturrock
#
#
# This file is part of LAMPSEQ.
# Copyright (C) 2025 Shane Sturrock
#
# LAMPSEQ is free software: you can redistribute it and/or modify
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


DIR=$(pwd)/$(dirname "$0")

rm -rf $DIR/01_fastq $DIR/02_qc $DIR/03_trim $DIR/04_align $DIR/05_results
