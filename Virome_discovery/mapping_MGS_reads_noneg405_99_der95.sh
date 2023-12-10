#!/bin/bash
#SBATCH --job-name=mapping_reads
#SBATCH --error=./err/mapping_reads.err
#SBATCH --output=./out/mapping_reads.out
#SBATCH --mem=8gb
#SBATCH --time=11:29:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

if [ -f ../Total_MGS/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq.gz ]; then
	gunzip ../Total_MGS/${SAMPLE_ID}_kneaddata_cleaned_pair_*.fastq.gz
fi

# PURGING ENVIRUMENT 
module purge 

# --- LOAD MODULES --- 
module load Python
module load Bowtie2
module load SAMtools

bowtie2 \
        -x ../clean_viral/dbs/viral_decontaminated/viral_decontaminated_db \
        -1 ../Total_MGS/${SAMPLE_ID}_kneaddata_cleaned_pair_1.fastq \
        -2 ../Total_MGS/${SAMPLE_ID}_kneaddata_cleaned_pair_2.fastq \
        -S ../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sam \
        -p ${SLURM_CPUS_PER_TASK} \
        --no-unal \
        --end-to-end \
        &> ../clean_viral/MGS_to_viral_decontaminated/alignment_log/${SAMPLE_ID}.bowtie2.log

samtools view \
	-S ../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sam \
	-b \
	-o ../Total_MGS/${SAMPLE_ID}_all_vir_alignments.bam

samtools sort \
	../Total_MGS/${SAMPLE_ID}_all_vir_alignments.bam \
	-o ../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sorted.bam

samtools index \
	../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sorted.bam

samtools mpileup \
	../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sorted.bam \
	-A \
	-f ../clean_viral/dbs/viral_decontaminated/viral_noneg405_99_der95_decontaminated.fasta \
	-o ../clean_viral/MGS_to_viral_decontaminated/breadth_cov/${SAMPLE_ID}.A.cov.txt

samtools idxstats \
	../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sorted.bam \
	> ../clean_viral/MGS_to_viral_decontaminated/counts_table/${SAMPLE_ID}.counts.idx.txt

samtools flagstat \
	../Total_MGS/${SAMPLE_ID}_all_vir_alignments.sorted.bam \
	 > ../clean_viral/MGS_to_viral_decontaminated/counts_table/${SAMPLE_ID}.flagstat.txt

rm ../Total_MGS/${SAMPLE_ID}_all_vir_alignments*

gzip ../Total_MGS/${SAMPLE_ID}_kneaddata_cleaned_pair_*.fastq

# PURGING ENVIRUMENT
module purge
