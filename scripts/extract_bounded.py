#!/usr/bin/env python3

# WildSEQ: a qPCR Sequencing Bioinformatics Pipeline
# Author: S. Sturrock
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



import sys, gzip, argparse
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from Bio.Seq import Seq

# Initialize parser
parser = argparse.ArgumentParser()

# Adding arguments
parser.add_argument("-f", "--flanks", help = "Flanks file name")
parser.add_argument("-i", "--input", help = "Input fastq file name")
parser.add_argument("-o", "--output", help = "Output fasta file name")
parser.add_argument("-v", "--verbose", action="store_true", help = "Verbose output")

# Read arguments from command line
args = parser.parse_args()

# Set Parameters

flank_set = False
input_set = False
output_set = False
verbose_set = args.verbose

# File containing start end stop sequences representing flank of
# target region
if args.flanks:
  flanking_seqs_file_name = args.flanks
  flank_set = True
# Input file
if args.input:
  input_file_name = args.input
  input_set = True
# Output file
if args.output:
  output_file_name = args.output
  output_set = True

if not flank_set:
  print("You need to provide a file containing the flanking pairs")
  exit()

if not input_set:
  print("You need to provide a file of reads to extract from")
  exit()

if not output_set:
  print("You need to provide an output file name to write to")
  exit()

encoding = guess_type(flanking_seqs_file_name)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

# Flanking seqs stored as matched pairs
flanking_seq_list = dict()
records = list(SeqIO.parse(flanking_seqs_file_name, "fasta"))
# Read through records in steps of two
for loop in range(0, len(records), 2):
    # Load the pairs of records
    flanking_seq1 = records[loop]
    flanking_seq2 = records[loop+1]
    flank1 = flanking_seq1.seq
    name1 = flanking_seq1.name
    parts1 = flanking_seq1.description.split(" ")
    len1 = int(parts1[1])
    flank2 = flanking_seq2.seq
    name2 = flanking_seq2.name
    parts2 = flanking_seq2.description.split(" ")
    len2 = int(parts2[1])
    if len1 != len2:
        print("Lengths don't match for",name1,"and",name2)
        print(len1,"!=",len2)
        exit(1)

    flanking_seq_list[name1, len1, flank1, name2, len2, flank2] = 1

encoding = guess_type(input_file_name)[1]  # uses file extension
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

# Find flanking regions in both strands

total_reads = 0

output_file = open(output_file_name, "w")

total = 0
with _open(input_file_name) as input_file:
    # read each sequence
    if verbose_set:
        print("\rTotal_reads:", total_reads, end = "")
    for record in SeqIO.parse(input_file, 'fastq'):
        total_reads += 1
        if total_reads % 1000 == 0 and verbose_set:
            print("\rTotal_reads:", total_reads, end = "")

        for flanking_name1, expected_len1, flanking_seq1, flanking_name2, expected_len2, flanking_seq2 in list(flanking_seq_list.keys()):
            # Find the target seq in the reverse complement sequence
            fs1loc_fwd = -1
            fs2loc_fwd = -1
            fs1loc_rev = -1
            fs2loc_rev = -1

            # Check forward direction
            if flanking_seq1 in record.seq:
                fs1loc_fwd = record.seq.find(flanking_seq1)
            if flanking_seq2 in record.seq:
                fs2loc_fwd = record.seq.find(flanking_seq2)

            # Print out the forward region as a read
            if (fs1loc_fwd != -1) and (fs2loc_fwd != -1):
                region = record.seq[fs1loc_fwd+len(flanking_seq1):fs2loc_fwd]
                if (len(region) > expected_len1-20) and (len(region) < expected_len1+20):
                  print(">Read_",total,"\n",region,sep="",file=output_file)
                  total += 1

            # Check reverse direction
            if flanking_seq1 in record.seq.reverse_complement():
                fs1loc_rev = record.seq.reverse_complement().find(flanking_seq1)
            if flanking_seq2 in record.seq.reverse_complement():
                fs2loc_rev = record.seq.reverse_complement().find(flanking_seq2)

            # Print out reverse region as a read
            if (fs1loc_rev != -1) and (fs2loc_rev != -1):
                region = record.seq.reverse_complement()[fs1loc_rev+len(flanking_seq1):fs2loc_rev]
                if (len(region) > expected_len1-20) and (len(region) < expected_len1+20):
                  print(">Read_",total,"\n",region,sep="",file=output_file)
                  total += 1

output_file.close()
if verbose_set:
    print("\rTotal_reads:", total_reads)
