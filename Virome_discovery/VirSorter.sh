#!/bin/bash
#SBATCH --job-name=VirSorter
#SBATCH --output=VirSorter.out
#SBATCH --error=VirSorter.err
#SBATCH --mem=4gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

module load Biopython
module load BioPerl
module load OpenBLAS
module load Anaconda3

source activate virsorter

wrapper_phage_contigs_sorter_iPlant.pl \
	-f ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta \
	--db 2 \
	--virome \
	--wdir ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/VirSorter \
	--ncpu 4 \
	--data-dir /home/umcg-sgarmaeva/virsorter-data
