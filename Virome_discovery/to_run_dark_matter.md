# Example of directory preparation to run dark_matter.sh
Since this criteria requires a lot of time even when running individual samples, I decided to run it in parallel by splitting the concatenated scaffolds in similar length fragments.

## Concatenating scaffolds per batch
```
mkdir /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/BATCH1/dark_matter
cd /scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/BATCH1/dark_matter
for i in $(cat ../../samples.batch1); do cat ../${i}/assembly_metaspades/${i}_scaffolds.min1000.fasta >> all_scaffolds.fasta; done
```
## Filtering by size
```
/data/umcg-sgarmaeva/scripts/filter_contigs.pl 3000 all_scaffolds.fasta > all_scaffolds.3k.fasta
rm all_scaffolds.fasta
```
## Calculating the number of nucleotides
```
module load QUAST
quast -o quast_out all_scaffolds.3k.fasta
```
Number of nucleotides is 1,066,550,882 (from quast report), and every split file should contain around 10,000,000 bases to finish running in the reasonable time
Thus, we need to split scaffolds to 107 fasta.

## Splitting 
```
/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/scripts/fastasplitn all_scaffolds.3k.fasta 107
```
[fastasplitn](https://github.com/ISUgenomics/common_scripts/blob/master/fastasplitn.c) is a modified and compiled script. Modifications: line #164:
```
#(void) strcpy(thetemplate,"frag%.3d");
(void) strcpy(thetemplate,"frag%.3d.fa");
```

ls *.fa | awk -F '.fa' '{print $1}' > fragment.list
