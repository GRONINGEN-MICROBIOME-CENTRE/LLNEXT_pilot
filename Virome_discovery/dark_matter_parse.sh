#!/bin/bash
#SBATCH --job-name=parsing
#SBATCH --output=parsing.out
#SBATCH --error=parsing.err
#SBATCH --partition=short
#SBATCH --mem=16gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1

BATCH_ID=$1
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

# PURGING ENVIRUMENT 
module purge

# PARSING IDs AND SOURCE

# Dark Matter:

cat ../BATCH${BATCH_ID}/dark_matter/frag*.txt > ../BATCH${BATCH_ID}/dark_matter/all_hits.txt

bash /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/nt_hits_filtering.sh \
	../BATCH${BATCH_ID}/dark_matter/all_hits.txt \
	90 \
	0
awk -F '\t' '$3 > 100' ../BATCH${BATCH_ID}/dark_matter/contigs_blastall_filtered.txt | \
	awk -F '\t' '{print $1}' | \
	sort | \
	uniq \
	> ../BATCH${BATCH_ID}/dark_matter/contigs_discard

awk -F '\t' '$3 > 90 && $4 > 100' ../BATCH${BATCH_ID}/dark_matter/all_hits.txt | \
	awk -F '\t' '{print $1}' | \
	sort | \
	uniq \
	> ../BATCH${BATCH_ID}/dark_matter/contigs_discard_add

cat ../BATCH${BATCH_ID}/dark_matter/contigs_discard ../BATCH${BATCH_ID}/dark_matter/contigs_discard_add | \
	sort | \
	uniq \
	> ../BATCH${BATCH_ID}/dark_matter/contigs_discard_all
 
grep '>' ../BATCH${BATCH_ID}/dark_matter/all_scaffolds.3k.fasta | \
	sed 's/>//g' | \
	sort \
	> ../BATCH${BATCH_ID}/dark_matter/all_contigs_3k

comm -23 ../BATCH${BATCH_ID}/dark_matter/all_contigs_3k ../BATCH${BATCH_ID}/dark_matter/contigs_discard_all > ../BATCH${BATCH_ID}/dark_matter/dark_matter_all

for SAMPLE_ID in $(awk -F '\t' '{print $2}' ../samples.batch${BATCH_ID}.renaming); do 

	grep ${SAMPLE_ID} ../BATCH${BATCH_ID}/dark_matter/dark_matter_all | \
	awk '{print $0 "\tdark_matter"}' \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/dark_matter_tidy 
done

