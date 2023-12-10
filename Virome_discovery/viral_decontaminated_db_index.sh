#!/bin/bash
#SBATCH --job-name=viral_decontaminated_db
#SBATCH --error=./err/viral_decontaminated_db.err
#SBATCH --output=./out/viral_decontaminated_db.out
#SBATCH --mem=6gb
#SBATCH --time=02:29:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

mkdir -p ../clean_viral/dbs/viral_decontaminated/
# PURGING ENV
module purge

# LOADING MODULES
module load Python
module load Bowtie2
module load SAMtools

# DB indexing
bowtie2-build \
	../clean_viral/scaffolds/viral_noneg405_99_der95_decontaminated.fasta \
	../clean_viral/dbs/viral_decontaminated/viral_decontaminated_db

samtools faidx \
	../clean_viral/scaffolds/viral_noneg405_99_der95_decontaminated.fasta

mv ../clean_viral/scaffolds/viral_noneg405_99_der95_decontaminated.fasta.fai ../clean_viral/dbs/viral_decontaminated/

awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' \
	../clean_viral/dbs/viral_decontaminated/viral_noneg405_99_der95_decontaminated.fasta.fai \
	> ../clean_viral/dbs/viral_decontaminated/viral_noneg405_99_der95_decontaminated.fasta.bed

mkdir -p ../clean_viral/VLP_to_viral_decontaminated/alignment_log
mkdir -p ../clean_viral/VLP_to_viral_decontaminated/breadth_cov
mkdir -p ../clean_viral/VLP_to_viral_decontaminated/counts_table

cp ../clean_viral/scaffolds/viral_noneg405_99_der95_decontaminated.fasta ../clean_viral/dbs/viral_decontaminated/


# PURGING ENV
module purge
