#!/bin/bash
#SBATCH --job-name=dark_matter
#SBATCH --output=dark_matter.out
#SBATCH --error=dark_matter.err
#SBATCH --mem=6gb
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

# PURGING ENVIRUMENT 
module purge

# local alignment against NT database
blastn \
	-query ../BATCH${BATCH_ID}/dark_matter/${SAMPLE_ID}.fa \
	-db /scratch/umcg-sgarmaeva/databases/nt/nt \
	-evalue 1e-10 \
	-outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
	-out ../BATCH${BATCH_ID}/dark_matter/${SAMPLE_ID}.txt \
	-num_threads 1


