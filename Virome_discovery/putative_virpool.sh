#!/bin/bash
#SBATCH --job-name=virpool
#SBATCH --output=virpool.out
#SBATCH --error=virpool.err
#SBATCH --partition=short
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

# PURGING ENVIRUMENT 
module purge
# Loading modules
module load seqtk

# Combining contigs IDs for all putative viruses
cat ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/*_tidy | \
	awk -F '\t' '{print $1}' | \
	sort | \
	uniq \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids

cat ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/*_tidy > ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_table_of_origin
cp ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_table_of_origin ../clean_viral/table_of_origin

cat ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/pVOG_stat >> ../clean_viral/table_of_origin/pVOG_stat_all

seqtk \
	subseq \
	-l60 \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/all_predicted_viral_ids \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta

module purge

cp ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/${SAMPLE_ID}_all_predicted_viral.fasta ../clean_viral/scaffolds
