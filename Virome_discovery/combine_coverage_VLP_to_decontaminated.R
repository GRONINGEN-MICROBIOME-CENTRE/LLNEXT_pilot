library(tidyverse)

temp = list.files(path="/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/VLP_to_viral_decontaminated/breadth_cov",pattern="*.A.cov.p.contig.txt",full.names=T)

myfiles = lapply(temp, 
		FUN = function(x){
		read.table(x, sep='\t', header=T)
		})

coverage_table <- myfiles %>% reduce(full_join, by='V1')

coverage_table[is.na(coverage_table)] <- 0

write.table(coverage_table, '../clean_viral/VLP_to_viral_decontaminated/coverage_table.txt', sep='\t', row.names=F, quote=F)
