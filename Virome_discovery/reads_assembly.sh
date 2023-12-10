#!/bin/bash
#SBATCH --job-name=rQC
#SBATCH --error=./err/rAs.err
#SBATCH --output=./out/rAs.out
#SBATCH --mem=48gb
#SBATCH --time=07:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

# --- MAKE FOLDERS ----
mkdir ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades

# --- WORKING IN $TMPDIR ---
mkdir -p ${TMPDIR}/${SAMPLE_ID}/clean_reads/
cp ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_*.fastq ${TMPDIR}/${SAMPLE_ID}/clean_reads/

# --- LOAD MODULES ---
module load Miniconda3/4.7.10
source activate /data/umcg-tifn/rgacesa/conda_dag3_assemblers

# --- METAGENOME ASSEMBLY ---
/data/umcg-tifn/rgacesa/conda_dag3_assemblers/bin/metaspades.py \
	-1 ${TMPDIR}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq \
	-2 ${TMPDIR}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq \
	-o ${TMPDIR}/${SAMPLE_ID}/assembly_metaspades \
	-m $((${SLURM_MEM_PER_NODE} / 1024)) \
	-t ${SLURM_CPUS_PER_TASK}

cp ${TMPDIR}/${SAMPLE_ID}/assembly_metaspades/contigs.fasta ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_metaspades_contigs.fa
cp ${TMPDIR}/${SAMPLE_ID}/assembly_metaspades/scaffolds.fasta ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_metaspades_scaffolds.fa
cp ${TMPDIR}/${SAMPLE_ID}/assembly_metaspades/spades.log ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_metaspades_log.log

source deactivate
