#!/bin/bash
#SBATCH --job-name=parsing
#SBATCH --output=parsing.out
#SBATCH --error=parsing.err
#SBATCH --partition=short
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy

# PURGING ENVIRUMENT 
module purge

# PARSING IDs AND SOURCE

# CrAss:
bash /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/blast_results_filtering.sh \
        ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss_UPD/min_1k_crAss.outfmt6.txt \
        50 \
        0.9 
awk -F '\t' '{print $1}' ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss_UPD/contigs_blastall_filtered.txt | \
        sort | uniq | awk '{print $0 "\tCrAss"}' > ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/crAss_UPD_tidy



