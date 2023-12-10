#!/bin/bash
#SBATCH --job-name=job3
#SBATCH --output=job3_%A.out
#SBATCH --mem=4gb
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




module purge; module load R; module list

Rscript build_VC_table.R

chmod 440 pilot_crAss_ProkaryoticViralRefSeq211.VC.txt
