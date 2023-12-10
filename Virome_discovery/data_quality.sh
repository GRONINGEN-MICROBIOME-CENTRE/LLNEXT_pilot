#!/bin/bash
#SBATCH --job-name=assembly_quality
#SBATCH --error=assembly_quality.err
#SBATCH --output=assembly_quality.out
#SBATCH --partition=short
#SBATCH --mem=8gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

# PURGING ENVIRUMENT 
module purge 


# CONTIGS TRIMMING AND RENAMING
/data/umcg-sgarmaeva/scripts/filter_contigs.pl 1000 \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_metaspades_scaffolds.fa > \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta

sed -i 's/>NODE/>'${SAMPLE_ID}'_NODE/g' ../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_scaffolds.min1000.fasta

# PURGING ENVIRUMENT 
module purge
# --- LOAD MODULES --- 
module load Python/3.8.2-GCCcore-9.3.0
module load QUAST

## assembly quality assessment
quast \
	-o ../BATCH${BATCH_ID}/${SAMPLE_ID}/quast_out \
	--min-contig 1000 \
        --threads ${SLURM_CPUS_PER_TASK} \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/assembly_metaspades/${SAMPLE_ID}_metaspades_scaffolds.fa

# PURGING ENVIRUMENT
module purge

gunzip ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/*.gz

# --- LOAD MODULES --- 
module load Python
module load Bowtie2

bowtie2 \
        -x /data/umcg-sgarmaeva/databases/cpn60db/cpn60db_apr_22/cpn60db \
        -1 ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq \
        -2 ../BATCH${BATCH_ID}/${SAMPLE_ID}/clean_reads/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq \
        -S ../BATCH${BATCH_ID}/${SAMPLE_ID}/${SAMPLE_ID}_cpn60db.sam \
        -p ${SLURM_CPUS_PER_TASK} \
        --no-unal \
        --end-to-end \
        &> ../BATCH${BATCH_ID}/${SAMPLE_ID}/bowtie2.cpn60db.log


module purge
