#!/bin/bash

#Description: Extract information about alignment rate to the database of viral contigs
#Created by: Sana
#Year: 2019 

DIR=$1
SAMPLE_LIST=$2

echo -e "SID\tperc_reads_aligned" >> ${DIR}/alignment_stat.txt
for SAMPLE in $(cat ${SAMPLE_LIST});
do
  	echo `basename ${SAMPLE}` >> tmp
	less ${DIR}/alignment_log/${SAMPLE}.bowtie2.log | grep "overall alignment rate" | awk -F " " '{print $1}' >> tmp
	less tmp | paste -s >> tmp2
        paste -d "\n" tmp2 >> ${DIR}/alignment_stat.txt
        rm tmp*
done
