#!/bin/bash
#SBATCH --job-name=viral_references
#SBATCH --output=viral_references.out
#SBATCH --error=viral_references.err
#SBATCH --partition=short
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/viral_refseq
mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss

# PURGING ENVIRUMENT 
module purge

# local alignment against Viral RefSeq database
blastn \
	-query ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta \
	-db /data/umcg-sgarmaeva/databases/viral_refseq/viral_refseq_mar_22/viral_refseq \
	-evalue 1e-10 \
	-outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
	-out ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/viral_refseq/min_1k_viral_refseq.outfmt6.txt \
	-num_threads 1


