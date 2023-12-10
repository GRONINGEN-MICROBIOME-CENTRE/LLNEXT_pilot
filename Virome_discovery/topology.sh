#!/bin/bash
#SBATCH --job-name=Topology
#SBATCH --output=Topology.out
#SBATCH --error=Topology.err
#SBATCH --mem=4gb
#SBATCH --time=01:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/topology_classificator/lastz-distrib-1.04.15
export PATH=$PATH:/data/umcg-sgarmaeva/tools/topology_classificator

mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/Topology

# PURGING ENVIRUMENT 
module purge

# LOADING MODULES
module load Biopython

# Topologically classify the sequences
find_circular.py \
	-i ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta \
	-m 1000 
mv /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/${SAMPLE_ID}_scaffolds.min1000.fasta_circular.fna ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/Topology
# PURGING ENVIRUMENT 
module purge

