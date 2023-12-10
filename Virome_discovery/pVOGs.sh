#!/bin/bash
#SBATCH --job-name=pVOGs
#SBATCH --output=pVOGs.out
#SBATCH --error=pVOGs.err
#SBATCH --mem=2gb
#SBATCH --time=02:29:00
#SBATCH --cpus-per-task=2

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs

# PURGING ENVIRUMENT 
module purge

# --- LOAD MODULES --- 
module load prodigal

# predicting ORFs
prodigal \
	-a ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.AA.fasta \
	-i ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta \
	-p meta \
	&> ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/prodigal.log

# PURGING ENVIRUMENT 
module purge

# --- LOAD MODULES --- 
module load HMMER

# pVOGs search
hmmsearch \
	-E 0.000001 \
	--tblout ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/${SAMPLE_ID}_scaffolds.min1000.AA.tblout \
	--cpu 4 \
	/data/umcg-sgarmaeva/pvogs/AllvogHMMprofiles/all_vogs.hmm \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.AA.fasta \
	&> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/hmmsearch.log
		
