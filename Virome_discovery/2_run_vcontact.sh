#!/bin/bash
#SBATCH --job-name=job2
#SBATCH --output=job2_%A.out
#SBATCH --mem=120gb
#SBATCH --partition=himem
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE
#SBATCH --get-user-env=L




# -------------------------------------- prepare input files -------------------------------------- #
prot_fa='pilot_crAss.AA.fasta'
cat pilot_viruses.AA.fasta > ${prot_fa}
cat crAss_phages.AA.fasta >> ${prot_fa}


g2g_csv='pilot_crAss.g2g.csv'
echo 'protein_id,contig_id,keywords' > ${g2g_csv}
grep '>' ${prot_fa} | perl -pe 's/^>(\S+)(_[0-9]+) .+$/\1\2,\1,None_provided/' >> ${g2g_csv}


chmod 440 pilot_crAss.*




# -------------------------------------- run vConTACT -------------------------------------- #
module purge; module load Anaconda3; module list
source activate vContact2; conda list

vcontact2 \
	--raw-proteins ${prot_fa} \
	--proteins-fp ${g2g_csv} \
    --db 'ProkaryoticViralRefSeq211-Merged' \
	--output-dir 'vcontact_output' \
	--c1-bin '/home/umcg-agulyaeva/SOFTWARE/cluster_one-1.0.jar' \
	--threads ${SLURM_CPUS_PER_TASK}

conda deactivate


chmod 750 vcontact_output
chmod 440 vcontact_output/*
