#!/bin/bash
#SBATCH --job-name=job1
#SBATCH --output=job1_%A.out
#SBATCH --mem=4gb
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# -------------------------------------- predict proteomes of pilot viruses -------------------------------------- #
pilot_fa=/data/umcg-llnext/tmp_for_Nastya/LLNEXT_pilot_fup_final/all_viral_noneg405_99_der95.fasta

md5sum ${pilot_fa}


module purge; module load prodigal; module list

prodigal \
    -p meta \
    -i ${pilot_fa} \
    -a pilot_viruses.AA.fasta


chmod 440 pilot_viruses.AA.fasta




# -------------------------------------- predict proteomes of crAss-like phages -------------------------------------- #
crAss_fa=/data/umcg-llnext/tmp_for_Nastya/LLNEXT_pilot_fup_final/sequences_all_crAss_50k_der99.fasta

md5sum ${crAss_fa}


prodigal \
    -p meta \
    -i ${crAss_fa} \
    -a crAss_phages.AA.fasta


chmod 440 crAss_phages.AA.fasta
