#!/bin/bash
#SBATCH --job-name=M4
#SBATCH --error=./err/M4.err
#SBATCH --output=./out/M4.out
#SBATCH --mem=20gb
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# --- LOAD MODULES ---
module purge 
module load Anaconda3/2022.05
conda activate /scratch/hb-tifn/condas/conda_biobakery4/

# --- RUN METAPHLAN4 --- 
mkdir -p ../BATCH5/${SAMPLE_ID}/metaphlan4

# prep node temp folders
mkdir -p ${TMPDIR}/${SAMPLE_ID}
mkdir -p ${TMPDIR}/${SAMPLE_ID}/metaphlan4
mkdir -p ${TMPDIR}/${SAMPLE_ID}/clean_reads

cat ../BATCH5/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq > ${TMPDIR}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_paired_merged.fastq
cat ../BATCH5/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq >> ${TMPDIR}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_paired_merged.fastq

metaphlan \
	${TMPDIR}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_paired_merged.fastq \
	--input_type fastq \
	--nproc ${SLURM_CPUS_PER_TASK} \
	-o ${TMPDIR}/${SAMPLE_ID}/metaphlan4/${SAMPLE_ID}_metaphlan.txt \
	--tmp_dir ${TMPDIR}/${SAMPLE_ID}/metaphlan4_tmp \
	--bowtie2db /scratch/hb-tifn/condas/conda_biobakery4/lib/python3.9/site-packages/metaphlan/metaphlan_databases/ \
	--force --unclassified_estimation 2>&1 | tee ../BATCH5/${SAMPLE_ID}/metaphlan4/${SAMPLE_ID}_metaphlan4.log

cp ${TMPDIR}/${SAMPLE_ID}/metaphlan4/* ../BATCH5/${SAMPLE_ID}/metaphlan4/
rm -r ${TMPDIR}/${SAMPLE_ID}

