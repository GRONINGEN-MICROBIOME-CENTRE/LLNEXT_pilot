#pVOGs_filtering_hmmsearch takes a cleaned output from hmmsearch of NR contigs against pVOG db 
#and returns a dataframe containing contig id validated as viral by pVOGs, its length, 
#number of pVOGs hist and number of pVOGs hits per 10kb length of the contig. 
#Current cut-offs are: Score >= 15, pVOG_hits >=3, pVOG_per_10kb >=2.
#In case of using HMMscan output, it is necessary to recombine columns in the initial
#table (so the input looks like 'Protein_ID', 'pVOG', 'Score', 'Contig_ID', 'Length').

args = commandArgs(trailingOnly=TRUE)

pVOGs_filtering_hmmsearch <- function(hmmsearch_output){
  colnames(hmmsearch_output)[1:5] <- c('Protein_ID', 'pVOG',  'Score', 'Contig_ID', 'Length')
  hmmsearch_output$pVOG <- as.character(hmmsearch_output$pVOG)
  hmmsearch_output$Protein_ID <- as.character(hmmsearch_output$Protein_ID)
  hmmsearch_output$Score <- as.numeric(hmmsearch_output$Score)
  hmmsearch_output$Contig_ID <- as.character(hmmsearch_output$Contig_ID)
  hmmsearch_output$Length <- as.numeric(hmmsearch_output$Length)
  hmmsearch_output <- hmmsearch_output[order(hmmsearch_output$Score, decreasing = T), ]
  hmmsearch_output <- hmmsearch_output[!duplicated(hmmsearch_output$Protein_ID), ]
  hmmsearch_output$genome_pVOG <- paste(hmmsearch_output$Contig_ID, hmmsearch_output$pVOG, sep = '_')
  hmmsearch_output <- hmmsearch_output[!duplicated(hmmsearch_output$genome_pVOG), ]
  hmmsearch_output$genome_pVOG <- NULL
  hmmsearch_output_tbl <- as.data.frame(table(hmmsearch_output$Contig_ID))
  colnames(hmmsearch_output_tbl)[1:2] <- c('Contig_ID', 'pVOG_hits')
  hmmsearch_output <- merge(hmmsearch_output, hmmsearch_output_tbl, by = 'Contig_ID', all.x = T)
  hmmsearch_output$pVOG_per_10kb <- (hmmsearch_output$pVOG_hits / hmmsearch_output$Length) * 10000
  pVOG_phage <- hmmsearch_output[which(hmmsearch_output$Length >= 1000 & hmmsearch_output$Score >= 15 & hmmsearch_output$pVOG_hits >= 3 & hmmsearch_output$pVOG_per_10kb >= 2), ]
  pVOG_phage <- pVOG_phage[order(pVOG_phage$Length, pVOG_phage$pVOG_hits, pVOG_phage$pVOG_per_10kb, decreasing = T), ]
  rownames(pVOG_phage) <- NULL
  pVOG_contigs <- pVOG_phage[!duplicated(pVOG_phage$Contig_ID), ]
  pVOG_contigs <- pVOG_contigs[,-c(2:4)]
  rownames(pVOG_contigs) <- NULL
  return(pVOG_contigs)
}
  
pVOG <- read.table(paste0(args[1], '/', 'pVOG_metadata.txt'), sep='\t', header = F)
pVOG_contigs <- pVOGs_filtering_hmmsearch(pVOG)
pVOG_contigs$Origin <- 'pVOGs'
write.table(pVOG_contigs[,c('Contig_ID', 'Origin')], paste0(args[1], '/','pVOG_contigs_tidy'), sep='\t', quote = F, row.names = F, col.names = F)
write.table(pVOG_contigs, paste0(args[1], '/','pVOG_stat'), sep='\t', quote = F, row.names = F, col.names = F)
