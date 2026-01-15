# WildSEQ: a qPCR Sequencing Bioinformatics Pipeline
# Authors: O. Yugovich
#
# This file is part of WildSEQ.
# Copyright (C) 2025 Olivia Yugovich
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




# DESCRIPTION:
# This script will read in the .bam files created by bowtie2, 
# Assign species ID based on position co-ordinates given in the species_config.csv file 
# Output .csv file with species calls and read counts


# DEPENDENCIES:
  # Ensure species_config.csv has the positions relevant to the reference sequences used
  # pos_start: starting position co-ordinates of species in reference sequence
  # pos_finish: finishing position co-ordinates of species in reference sequence
  # species_id: species name related to pos_start and pos_finish co-ordinates

# load packages - quiet
load_packages_quietly <- function(pkgs) {
  invisible(lapply(pkgs, function(pkg) {
      suppressPackageStartupMessages(
        suppressWarnings(library(pkg, character.only = TRUE, quietly = TRUE)))
    })
  )
}

pkgs <- c("tidyverse", "dplyr", "Rsamtools", "here")
load_packages_quietly(pkgs)

# FOR BAM IMPORT and SPECIES_CONFIG IMPORT:
  # get bash command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  # get named arguments for INDIR and SCRIPT_DIR:
  get_arg <- function(flag, default = NULL) {
    match <- which(args == flag)
   if (length(match) == 1 && match < length(args)) {
      return(args[match + 1])
   } else {
     return(default)
   }
  }

  # get INDIR and SCRIPT_DIR from command arguments
  inputdir <- normalizePath(get_arg("--indir"), mustWork = TRUE)
  scriptdir <- normalizePath(get_arg("--scriptdir"), mustWork = TRUE)


  # FOR SPECIES_CONFIG IMPORT:
  # this file detonates the position on the concatenated reference sequence related to species
  config_path <- file.path(scriptdir, "species_id", "species_config.csv")
  config <- read_csv(config_path, show_col_types = FALSE)


  # FOR BAM FILE IMPORTS:
  bam_path <- file.path(inputdir, "04_align", "bam")
  setwd(bam_path)
  sum.files = list.files(pattern="\\.bam$")
  bam.import <- file.path(bam_path, sum.files)
# REQUIRED FILES HAVE NOW BEEN IMPORTED - continue analysis 
  
# ASSIGN SPECIES ID BASED ON POSITION IN REF_SEQ.FASTA
# assign_species function: use species_config parameters to determine species calls
assign_species <- function(pos) {
  if (is.na(pos)) return(NA)
  
  for (i in 1:nrow(config)) {
    pos_start <- config$pos_start[i]
    pos_finish <- config$pos_finish[i]
    
    if (pos >= pos_start && pos <= pos_finish) {
      return(config$species_id[i])
    }
    
  }
  return(NA)
}

# define columns to extract from bams
param <- ScanBamParam(what = c("qname", "pos"))
process_bam_file <- function(bamfile) {
  aln <- scanBam(bamfile, param = param)[[1]]
  
  # Skip bam if file contains no alignments
  if (length(aln$qname) == 0 || length(aln$pos) == 0) {
    message("Skipping empty file: ", bamfile)
    return(NULL)
  }
  
  df <- data.frame(
    qname = aln$qname,
    pos = aln$pos,
    stringsAsFactors = FALSE
  )
  
  # Remove rows with NA in 'pos'
  df <- df[!is.na(df$pos), ]
  
  if (nrow(df) == 0) {
    message("No valid positions in file: ", bamfile)
    return(NULL)
  }
  
  # Assign species
  df$species_ID <- sapply(df$pos, assign_species)
  
  # Keep only rows with a species match
  df <- df[!is.na(df$species_ID), ]
  
  if (nrow(df) == 0) {
    message("No matching species regions in file: ", bamfile)
    return(NULL)
  }
  
  # Assign barcode name from bam file name
  df$barcode <- tools::file_path_sans_ext(basename(bamfile))
  
  return(df)
}


# Apply to all bam files
species_assigned_list <- lapply(bam.import, process_bam_file)

# Combine all processed bams into one full data frame
combined.bam.df <- do.call(rbind, species_assigned_list)

# count number of reads per species ID result and print in "reads" column
species_calls <- aggregate(qname ~ barcode + species_ID, data = combined.bam.df, FUN = length)
colnames(species_calls)[3] <- "reads"

# reorder species_calls by barcode then species_ID
species_calls <- species_calls %>% 
  arrange(barcode, species_ID) %>% 
  mutate(species_ID = ifelse(species_ID == "M.primigenius", "Mammoth", species_ID)) # change species ID to just Mammoth as cant be resolved

# name the df based on bash inputdir
sampleid <- basename(inputdir)
df_name <- paste0(sampleid, "_species_calls")

# assign the df name to species_calls
assign(df_name, species_calls)

#write df to csv located in indir folder
output_path <- file.path(inputdir, "05_results", paste0(df_name, ".csv"))
write_csv(species_calls, output_path)

# print confirmation message
cat("Species calls output written to csv, located in:", output_path)



##########
# font setup script for pdf rendering
# create font_setup.tex
write_font_setup <- function(scriptdir) {
  fonts_dir <- file.path(scriptdir, "species_id", "fonts")
  output_path <- file.path(scriptdir, "species_id", "font_setup.tex")
  
  font_setup_text <- sprintf(
    '\\usepackage{fontspec}
\\setmainfont{Nunito-Regular}[Path=%s/, Extension=.ttf]
\\setsansfont{Nunito-Regular}[Path=%s/, Extension=.ttf]
\\setmonofont{Nunito-Regular}[Path=%s/, Extension=.ttf]
\\newfontfamily\\boldfont{Nunito-Bold}[Path=%s/, Extension=.ttf]
\\renewcommand{\\familydefault}{\\sfdefault}',
    fonts_dir, fonts_dir, fonts_dir, fonts_dir
  )
  
  writeLines(font_setup_text, con = output_path)
}

write_font_setup(scriptdir)
