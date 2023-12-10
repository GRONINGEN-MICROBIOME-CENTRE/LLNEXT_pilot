#!/bin/bash
#SBATCH --job-name=viral_references
#SBATCH --output=viral_references.out
#SBATCH --error=viral_references.err
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss_UPD

# PURGING ENVIRUMENT 
module purge

# local alignment against new crAss database

blastn \
        -query ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta \
        -db /data/umcg-sgarmaeva/databases/all_crAssvirales_Oct22/sequences_all_crAss_50k_der99_db \
        -evalue 1e-10 \
        -outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle' \
        -out ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss_UPD/min_1k_crAss.outfmt6.txt \
        -num_threads 1

