#!/bin/bash
#SBATCH --job-name=parsing
#SBATCH --output=parsing.out
#SBATCH --error=parsing.err
#SBATCH --partition=short
#SBATCH --mem=4gb
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1

SAMPLE_ID=$1
echo "SAMPLE_ID=${SAMPLE_ID}"

BATCH_ID=$2
echo "BATCH_ID=${BATCH_ID}"

export PATH=$PATH:/data/umcg-sgarmaeva/tools/ncbi-blast-2.13.0+/bin

mkdir -p ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy

# PURGING ENVIRUMENT 
module purge

# PARSING IDs AND SOURCE

# VirSorter:
grep -v "#" ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/VirSorter/VIRSorter_global-phage-signal.csv | \
	cut -f1 -d "," | \
	sed s/VIRSorter_// | \
	cut -f 1 -d "-" | \
	sed s/_/./10 | \
	sort | uniq | awk '{print $0 "\tVirSorter"}' > ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/virsorter_tidy

# pVOGs:
sed -e '1,3d' < ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/${SAMPLE_ID}_scaffolds.min1000.AA.tblout | \
	 head -n -10 | awk -F ' ' '{print $1"\t"$3"\t"$6}' \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp
awk -F '\t' '{print $1}' ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp | \
	awk 'BEGIN {FS="_"; OFS="_"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp1
awk -F '_' '{print $8}' ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp1 \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp3
paste ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp1 \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp3 \
	> ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/pVOG_metadata.txt
rm ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/tmp*

module load R
Rscript /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/pVOG_filtering_function.R ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/

# PURGING ENVIRUMENT 
module purge
mv ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/pVOGs/pVOG_contigs_tidy ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/

# Viral RefSeq:
bash /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/blast_results_filtering.sh \
	../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/viral_refseq/min_1k_viral_refseq.outfmt6.txt \
	50 \
	0.9 
awk -F '\t' '{print $1}' ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/viral_refseq/contigs_blastall_filtered.txt | \
	sort | uniq | awk '{print $0 "\tVir_Ref"}' > ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/refseq_tidy

# CrAss:
bash /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/blast_results_filtering.sh \
        ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss/min_1k_crAss.outfmt6.txt \
        50 \
        0.9 
awk -F '\t' '{print $1}' ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/crAss/contigs_blastall_filtered.txt | \
        sort | uniq | awk '{print $0 "\tCrAss"}' > ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/crAss_tidy

# Topology:
grep -E '>' ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/Topology/${SAMPLE_ID}_scaffolds.min1000.fasta_circular.fna | \
	awk -F '>' '{print $2}' | \
	awk '{print $0 "\tTopology"}' > ../BATCH${BATCH_ID}/${SAMPLE_ID}/virome_discovery/tidy/topology_tidy 


