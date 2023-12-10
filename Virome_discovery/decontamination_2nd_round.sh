#!/bin/bash
#SBATCH --job-name=decontamination
#SBATCH --error=./err/decontamination.err
#SBATCH --output=./out/decontamination.out
#SBATCH --partition=short
#SBATCH --mem=3gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1
#SBATCH --open-mode=truncate

# PURGING ENV
module purge

# LOADING MODULES
module load R

Rscript table_of_origin.R

mkdir /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination
cp ../clean_viral/table_of_origin/table_of_origin_noneg405_99_der95 ../clean_viral/decontamination
cp ../agulyaeva_LLNEXTpilot_vConTACT/pilot_crAss_ProkaryoticViralRefSeq211.VC.txt ../clean_viral/decontamination
cp ../clean_viral/ribosomal_proteins/noneg405_99_der95/RbS_proteins_only_50_90 ../clean_viral/decontamination
cp ../clean_viral/rrna_PILOT_vir/PILOT_vir_rrna_contigs.ids ../clean_viral/decontamination

Rscript decontamination_rbs_rrna.R

# PURGING ENV
module purge

# LOADING MODULES
module load seqtk

seqtk \
	subseq \
        -l60 \
	../clean_viral/scaffolds/all_viral_noneg405_99_der95.fasta \
	../clean_viral/decontamination/viral_IDs_noneg405_99_der95_decontaminated \
	> ../clean_viral/scaffolds/viral_noneg405_99_der95_decontaminated.fasta

# PURGING ENV
module purge
