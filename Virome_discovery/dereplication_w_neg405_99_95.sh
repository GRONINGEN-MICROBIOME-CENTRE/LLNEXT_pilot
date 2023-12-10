#!/bin/bash
#SBATCH --job-name=dereplicating
#SBATCH --output=./out/dereplicating_w_neg405_99_95.out
#SBATCH --error=./err/dereplicatin_w_neg405_99_95.err
#SBATCH --mem=16gb
#SBATCH --time=08:29:00
#SBATCH --cpus-per-task=4


# PURGING ENVIRUMENT 
module purge

# LOADING MODULES
module load Anaconda3

conda activate checkv

cat ../clean_viral/scaffolds/LN_*_all_predicted_viral.fasta ../clean_viral/scaffolds/LN_7C08_VL_405_scaffolds.min1000.fasta > ../clean_viral/scaffolds/all_viral_w_neg405.fasta

makeblastdb \
	-in ../clean_viral/scaffolds/all_viral_w_neg405.fasta \
	-dbtype nucl \
	-out ../clean_viral/dereplication/all_viral_w_neg405_db

blastn \
	-query ../clean_viral/scaffolds/all_viral_w_neg405.fasta \
	-db ../clean_viral/dereplication/all_viral_w_neg405_db \
	-outfmt '6 std qlen slen' \
	-max_target_seqs 10000 \
	-out ../clean_viral/dereplication/all_viral_w_neg405_blast.tsv \
	-num_threads 4

module load Python
module load Biopython

python /data/umcg-sgarmaeva/tools/checkv-scripts/anicalc.py \
	-i ../clean_viral/dereplication/all_viral_w_neg405_blast.tsv \
	-o ../clean_viral/dereplication/all_viral_w_neg405_ani.tsv

python /data/umcg-sgarmaeva/tools/checkv-scripts/aniclust.py \
	--fna ../clean_viral/scaffolds/all_viral_w_neg405.fasta \
	--ani ../clean_viral/dereplication/all_viral_w_neg405_ani.tsv \
	--out ../clean_viral/dereplication/all_viral_w_neg405_clusters_99.tsv \
	--min_ani 99 \
	--min_tcov 85 \
	--min_qcov 0

grep \
	-vE "LN_7C08_VL_405" ../clean_viral/dereplication/all_viral_w_neg405_clusters_99.tsv | \
	awk -F '\t' '{print $2}' | \
	tr , '\n' | sort | uniq \
	> ../clean_viral/dereplication/all_viral_noneg405_99_IDs

module load seqtk

seqtk \
	subseq \
	-l60 \
	../clean_viral/scaffolds/all_viral_w_neg405.fasta \
	../clean_viral/dereplication/all_viral_noneg405_99_IDs \
	> ../clean_viral/scaffolds/all_viral_noneg405_99.fasta

makeblastdb \
       -in ../clean_viral/scaffolds/all_viral_noneg405_99.fasta \
       -dbtype nucl \
       -out ../clean_viral/dereplication/all_viral_noneg405_99_db

blastn \
       -query ../clean_viral/scaffolds/all_viral_noneg405_99.fasta \
       -db ../clean_viral/dereplication/all_viral_noneg405_99_db \
       -outfmt '6 std qlen slen' \
       -max_target_seqs 10000 \
       -out ../clean_viral/dereplication/all_viral_noneg405_99_blast.tsv \
       -num_threads 4 

module load Python
module load Biopython

python /data/umcg-sgarmaeva/tools/checkv-scripts/anicalc.py \
       -i ../clean_viral/dereplication/all_viral_noneg405_99_blast.tsv \
       -o ../clean_viral/dereplication/all_viral_noneg405_99_ani.tsv

python /data/umcg-sgarmaeva/tools/checkv-scripts/aniclust.py \
        --fna ../clean_viral/scaffolds/all_viral_noneg405_99.fasta \
        --ani ../clean_viral/dereplication/all_viral_noneg405_99_ani.tsv \
        --out ../clean_viral/dereplication/all_viral_noneg405_99_der95.tsv \
        --min_ani 95 \
        --min_tcov 85 \
        --min_qcov 0

awk -F '\t' '{print $1}' \
	../clean_viral/dereplication/all_viral_noneg405_99_der95.tsv \
	> ../clean_viral/dereplication/all_viral_noneg405_99_der95_IDs

seqtk \
        subseq \
        -l60 \
        ../clean_viral/scaffolds/all_viral_noneg405_99.fasta \
        ../clean_viral/dereplication/all_viral_noneg405_99_der95_IDs \
        > ../clean_viral/scaffolds/all_viral_noneg405_99_der95.fasta

