#!/bin/bash
#SBATCH --job-name=rQC
#SBATCH --error=./err/rQC.err
#SBATCH --output=./out/rQC.out
#SBATCH --mem=16gb
#SBATCH --time=07:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

# based on microbiome_profiling scripts from https://github.com/GRONINGEN-MICROBIOME-CENTRE/DMP

# --- MAKE FOLDERS ---- 
mkdir ../BATCH${BATCH_ID}/${SAMPLE_ID}
mkdir ../BATCH${BATCH_ID}/${SAMPLE_ID}/qc_preclean
mkdir ../BATCH${BATCH_ID}/${SAMPLE_ID}/qc_postclean
mkdir ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads

# --- LOAD MODULES ---
module load FastQC/0.11.7-Java-1.8.0_162
module load Miniconda3/4.7.10
source activate /data/umcg-tifn/rgacesa/conda_BB

# --- RAW READS QC ---
echo "Running FastQC on unclean reads"
fastqc -t ${SLURM_CPUS_PER_TASK} -q -o ../BATCH${BATCH_ID}/${SAMPLE_ID}/qc_preclean ../BATCH${BATCH_ID}/${SAMPLE_ID}_1.fq.gz
fastqc -t ${SLURM_CPUS_PER_TASK} -q -o ../BATCH${BATCH_ID}/${SAMPLE_ID}/qc_preclean ../BATCH${BATCH_ID}/${SAMPLE_ID}_2.fq.gz

# --- WORKING IN $TMPDIR ---
mkdir -p ${TMPDIR}/${SAMPLE_ID}/filtering_data/
cp ../BATCH${BATCH_ID}/${SAMPLE_ID}_1.fq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/
cp ../BATCH${BATCH_ID}/${SAMPLE_ID}_2.fq.gz ${TMPDIR}/${SAMPLE_ID}/filtering_data/

# --- RUN ADAPTER TRIMMING ---- 
echo "Running BBDUK" 
bbduk.sh \
        in1=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_1.fq.gz \
        in2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_2.fq.gz \
        out1=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_bbdukout_1.fq.gz \
        out2=${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_bbdukout_2.fq.gz \
        ref=/data/umcg-llnext/MGS_analysis/adapters/next_cohort_adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo 2>&1 \
        threads=${SLURM_CPUS_PER_TASK} \
        -Xmx$((${SLURM_MEM_PER_NODE} / 1024))g | tee -a ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_bbduk.log

source deactivate

# --- RUN KNEADDATA ---- 
echo "Running Kneaddata" 
source activate /data/umcg-tifn/rgacesa/conda_dag3_v3

kneaddata \
	--input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_bbdukout_1.fq.gz \
	--input ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_bbdukout_2.fq.gz \
	--threads ${SLURM_CPUS_PER_TASK} \
	--processes 4 \
	--output-prefix ${SAMPLE_ID}_kneaddata \
	--output ${TMPDIR}/${SAMPLE_ID}/filtering_data \
	--log ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_kneaddata.log \
	-db /data/umcg-tifn/rgacesa/dag3_pipeline_v3_dbs/GRCh38p13 \
	--trimmomatic /data/umcg-tifn/rgacesa/conda_dag3_v3/share/trimmomatic-0.39-1/trimmomatic.jar \
        --run-trim-repetitive \
        --fastqc fastqc \
        --sequencer-source none \
        --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" \
        --bypass-trf

cp ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_1.fastq ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq
cp ${TMPDIR}/${SAMPLE_ID}/filtering_data/${SAMPLE_ID}_kneaddata_paired_2.fastq ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq

rm -r ${TMPDIR}/${SAMPLE_ID}/filtering_data/

if [ -f ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq ]; then
	rm ../BATCH${BATCH_ID}/${SAMPLE_ID}_1.fq.gz
	rm ../BATCH${BATCH_ID}/${SAMPLE_ID}_2.fq.gz
	
	# --- CLEAN READS QC --- 
	echo "Running FastQC on cleaned reads"
	fastqc -t ${SLURM_CPUS_PER_TASK} -q -o ../BATCH${BATCH_ID}/${SAMPLE_ID}/qc_postclean ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq
	fastqc -t ${SLURM_CPUS_PER_TASK} -q -o ../BATCH${BATCH_ID}/${SAMPLE_ID}/qc_postclean ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq
	
	# --- LAUNCHING ASSEMBLY ---
	sbatch --output ./out/${SAMPLE}_rAs.out --error ./err/${SAMPLE}_rAs.err --job-name ${SAMPLE} reads_assembly.sh ${SAMPLE} ${BATCH_ID}
fi

