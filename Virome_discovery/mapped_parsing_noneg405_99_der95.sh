#!/bin/bash
#SBATCH --job-name=mapping_reads_parsing
#SBATCH --error=./err/mapping_reads_parsing.err
#SBATCH --output=./out/mapping_reads_parsing.out
#SBATCH --mem=16gb
#SBATCH --time=01:29:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# PURGING ENVIRUMENT 
module purge 

awk \
	-F '\t' '{print $1"\t"$4}' \
	../clean_viral/VLP_to_viral_decontaminated/breadth_cov/${SAMPLE_ID}.A.cov.txt \
	> ../clean_viral/VLP_to_viral_decontaminated/breadth_cov/${SAMPLE_ID}.A.cov.p.base.txt

# LOADING MODULES
module load R

Rscript breadth_coverage.R ../clean_viral/VLP_to_viral_decontaminated/breadth_cov/${SAMPLE_ID}

# PURGING ENVIRUMENT
module purge
