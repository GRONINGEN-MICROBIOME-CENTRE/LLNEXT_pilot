#!/bin/bash
#SBATCH --job-name=cogs_search
#SBATCH --output=cogs_search.out
#SBATCH --error=cogs_search.err
#SBATCH --mem=4gb
#SBATCH --time=6:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

# PURGING ENVIRUMENT 
module purge

# LOADING MODULES
module load prodigal

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

# predicting ORFs
prodigal \
        -a ../clean_viral/ribosomal_proteins/noneg405_99_der95/${SAMPLE_ID}.faa \
        -i ../clean_viral/ribosomal_proteins/noneg405_99_der95/${SAMPLE_ID}.fa \
        -p meta \
        &> ../clean_viral/ribosomal_proteins/noneg405_99_der95/prodigal_log/${SAMPLE_ID}_prodigal.log

# local alignment against a subset of COGs database
blastp \
	-query ../clean_viral/ribosomal_proteins/noneg405_99_der95/${SAMPLE_ID}.faa \
 	-db /data/umcg-sgarmaeva/databases/COG2020/COG2020_RbS/cogs_rbs_db \
	-evalue 1e-10 \
	-outfmt '6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle qcovs' \
	-out ../clean_viral/ribosomal_proteins/noneg405_99_der95/${SAMPLE_ID}.RbS.outfmt6.txt \
	-num_threads 1

