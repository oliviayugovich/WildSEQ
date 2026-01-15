#!/usr/bin/env bash

# WildSEQ: a qPCR Sequencing Bioinformatics Pipeline
# Authors: O. Yugovich, S. Sturrock

# This file is part of WildSEQ.
# Copyright (C) 2025 Olivia Yugovich and Shane Sturrock
#
# WildSEQ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 

# Fail on errors so the script doesn't continue
#set -e

# Set programme arguments and usage of programme
usage() 
{
echo "Usage: 
   --indir <value> set input directory (defaults to current working directory) #wd, --input-path
   --threads <value> set number of threads to use (defaults to 8)
   --dryrun Set dryrun flag for testing without running (default=false) # dryrun for testing
   --debug Set debug flag for verbose output (default=false unless dryrun set) # show errors
   --help Print this message # print help
  
  Note: Script to run WildSEQ analysis pipeline"
exit 2
}

# Call getopt to validate the provided input, set argument options
OPTS=$(getopt -o i:,t:,h --long indir:,threads:,dryrun,debug,help -- "$@")
VALID_ARGS=$?

# Make new lines the only separator, if no valid argument is given, print usage
IFS=$'\n'
if [ "$VALID_ARGS" != "0" ]; then 
  usage
fi
#echo "$OPTS"

# Define the variables for options INDIR, DRYRUN and DEBUG for programme
INDIR=`pwd`
DRYRUN=false
DEBUG=false
#THREADS=8
eval set -- "$OPTS"
while :
do
  case "$1" in
    --indir ) INDIR=$(echo "$2" | sed 's:/*$::'); shift 2 ;;
#    --threads ) THREADS="$2"; shift 2 ;;
    --dryrun ) DRYRUN=true; DEBUG=true; shift ;;
    --debug ) DEBUG=true; shift ;;
    -h | --help )
       usage
       ;;
    -- ) shift; break ;;
    *) echo "Unexpected argument: $1 - this should not happen."
       usage
       ;;
  esac
done


if [[ "$DRYRUN" = true ]]; then
  echo "Dryrun mode, no code will be run"
fi

if [[ "$DEBUG" = true ]]; then
  echo "--indir $INDIR"
  echo "--dryrun $DRYRUN"
fi


# 01 - Analysis Set Up
  # Make all analysis subdirectories
echo "Setting up WildSEQ workflow for analysis"

mkdir -p ${INDIR}/01_fastq/
mkdir -p ${INDIR}/02_qc/fastqc/
mkdir -p ${INDIR}/02_qc/multiqc/
mkdir -p ${INDIR}/03_trim/
mkdir -p ${INDIR}/04_align/sam/
mkdir -p ${INDIR}/04_align/bam/
mkdir -p ${INDIR}/05_results/

# Check barcode naming, clean up if necessary
# Copy raw fastq in INDIR 00_raw/ to 01_fastq/
RAW_DIR=${INDIR}/00_raw_data
FASTQ_DIR=${INDIR}/01_fastq

# Loop over all fastq files in the raw data directory
for file in "$RAW_DIR"/*.fastq; do
    base=$(basename "$file")
    
    # Check if the filename is already in the correct format: barcodeNN.fastq
    if [[ "$base" =~ ^barcode([0-9]{2})\.fastq$ ]]; then
        cp "$file" "$FASTQ_DIR/$base"
	if [[ "$DEBUG" = true ]]; then
            echo "Copied: $base"
	fi
    else
        # Extract barcodeNN from anywhere in the filename
        barcode=$(echo "$base" | sed -n 's/.*\(barcode[0-9][0-9]\).*/\1.fastq/p')
        
        if [[ -n "$barcode" ]]; then
            cp "$file" "$FASTQ_DIR/$barcode"
	    if [[ "$DEBUG" = true ]]; then
                echo "Renamed and copied: $base -> $barcode"
	    fi
        else
	    if [[ "$DEBUG" = true ]]; then
                echo "Skipped (invalid file name): $base"
	    fi
        fi
    fi
done

echo "WildSEQ workflow set-up completed successfully"

# 02 - Quality Control
  # FASTQC each barcode*.fastq file separately
echo "Performing FASTQC for per barcode QC"
if [[ "$DEBUG" = true ]]; then
    #fastqc --threads $THREADS -o ${INDIR}/02_qc/fastqc/ ${INDIR}/01_fastq/*.fastq | grep -v "null"
    fastqc -o ${INDIR}/02_qc/fastqc/ --memory 4096 ${INDIR}/01_fastq/*.fastq | grep -v "null"
else 
    #fastqc --threads $THREADS -o ${INDIR}/02_qc/fastqc/ ${INDIR}/01_fastq/*.fastq 2>/dev/null | grep "Analysis"
    fastqc -o ${INDIR}/02_qc/fastqc/ --memory 4096 ${INDIR}/01_fastq/*.fastq 2>/dev/null | grep "Analysis"
fi
echo "FASTQC on fastq files complete"

# MULTIQC all fastq files to analyse whole run QC metrics

echo "Performing MULTIQC for run QC"
MULTIQC_FLAGS="--quiet"
if [[ "$DEBUG" = true ]]; then
    MULTIQC_FLAGS=""
fi
multiqc $MULTIQC_FLAGS --outdir ${INDIR}/02_qc/multiqc/ --filename ${INDIR}_multiqc.html ${INDIR}/02_qc/fastqc/

# copy multiqc report to 05_results
cp ${INDIR}/02_qc/multiqc/${INDIR}_multiqc.html ${INDIR}/05_results/
echo "MULTIQC Report (.html) for ${INDIR} ready for review in ${INDIR}/05_results"

# 03 - Trimming & Alignment
# Set filepath variables
SCRIPT_DIR=$(pwd)/$(dirname "$0")
INPUT_DIR=${INDIR}/01_fastq
OUTPUT_DIR=${INDIR}/03_trim
ALIGN_DIR=${INDIR}/04_align

# Search for flanking regions and trim 
echo "Starting flank search and trimming"
EXTRACT_FLAGS=""
if [[ "$DEBUG" = true ]]; then
    EXTRACT_FLAGS="--verbose"
fi
for file in ${INDIR}/01_fastq/*.fastq; do
    BASENAME=$(echo $(basename $file) | sed 's/[.]fastq//g')
    echo -ne \\r"Trimming ${BASENAME}"
    if [[ "$DEBUG" = true ]]; then
        echo ""
    fi 
    ${SCRIPT_DIR}/extract_bounded.py -f ${SCRIPT_DIR}/flanks.fasta -i ${INPUT_DIR}/${BASENAME}.fastq -o ${OUTPUT_DIR}/${BASENAME}.fasta $EXTRACT_FLAGS
done
echo
echo "Trimming completed successfully"

# 04 - Alignment
# align using bowtie2
BOWTIE2_FLAGS="--quiet"
if [[ "$DEBUG" = true ]]; then
    BOWTIE2_FLAGS=""
fi
echo "Starting alignments with bowtie2"
for file in ${OUTPUT_DIR}/*.fasta; do
    BASENAME=$(echo $(basename $file) | sed 's/[.]fasta//g')
    echo -ne \\r"Aligning ${BASENAME} to WildSEQ species"
    if [[ "$DEBUG" = true ]]; then
        echo ""
    fi 
    bowtie2-align-s -I 0 -X 800 -p 12 --very-fast -f -x ${SCRIPT_DIR}/ref_seqs/ref_seqs -U ${OUTPUT_DIR}/${BASENAME}.fasta -S ${ALIGN_DIR}/sam/${BASENAME}.sam $BOWTIE2_FLAGS
done
echo
echo "Alignments with bowtie2 completed successfully"

# Convert sam to bam using samtools
for file in ${ALIGN_DIR}/sam/*.sam; do                             
    BASENAME=$(echo $(basename $file) | sed 's/[.]sam//g')
    if [[ "$DEBUG" = true ]]; then
      echo "Converting $BASENAME from sam to bam"
    fi
    samtools view -S -b ${ALIGN_DIR}/sam/${BASENAME}.sam > ${ALIGN_DIR}/bam/${BASENAME}.bam
done
echo "Converting files to bam outputs complete"
echo "WildSEQ analysis pipeline completed successfully - ready for R analysis"


# 05 - Results
  # set file path variables
INDIR=$(realpath "$INDIR")
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RESULTS_DIR=${INDIR}/05_results/
INDIR_NAME=$(basename "$INDIR")
FONTS_DIR="${SCRIPT_DIR}/species_id/fonts"

# take alignments (bams) and output csv with species_ID information
Rscript ${SCRIPT_DIR}/species_id/species_id.R --indir ${INDIR} --scriptdir ${SCRIPT_DIR}

# create Rmarkdown report based on output csv
Rscript -e "rmarkdown::render(
  input = '${SCRIPT_DIR}/species_id/WildSEQ_run_summary.Rmd', 
  params = list(
  indir = '${INDIR}', 
  results_dir = '${RESULTS_DIR}',
  fonts_dir = '${FONTS_DIR}'),
  output_file = paste0('WildSEQ_', '${INDIR_NAME}', '_run_summary.pdf'),
  output_dir = '${RESULTS_DIR}')"

# remove temp font_setup.tex file
rm -f "${SCRIPT_DIR}/species_id/font_setup.tex"

echo "
Run completed, output files located in ${RESULTS_DIR}"


