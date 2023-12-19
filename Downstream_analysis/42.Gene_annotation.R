setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore gene content of selected bifidophages
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################

library(tidyverse)
library(gggenomes)

##############################
# Input data
##############################
prot_pos <- read.delim("02.CLEAN_DATA/PROTEIN_POS_viral_noneg405_99_der95_decontaminated.txt", 
                       header = FALSE , sep="#")
colnames(prot_pos) <- c("protein", "coord_start", "coord_stop", "coord_strand", "annotation_info")
prot_pos$protein <- gsub('>', '', prot_pos$protein)
prot_pos$protein <- gsub(' ', '', prot_pos$protein)

amg_annotation <- read_tsv('02.CLEAN_DATA/VIBRANT_annotations_viral_noneg405_99_der95_decontaminated.tsv')

L34922_LS1_gene_annotation <- merge(prot_pos[grep('LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224', prot_pos$protein),], 
                                    amg_annotation[amg_annotation$scaffold=='LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224',])
colnames(L34922_LS1_gene_annotation)[c(2,3,4,6)] <- c( 'start', 'end', 'strand', 'seq_id')

L37775_LS1_gene_annotation <- merge(prot_pos[grep('LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438', prot_pos$protein),], 
                                    amg_annotation[amg_annotation$scaffold=='LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438',])
colnames(L37775_LS1_gene_annotation)[c(2,3,4,6)] <- c( 'start', 'end', 'strand', 'seq_id')

#### prophage from Bifidum
prot_pos_pro <- read.delim("03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/PROTEIN_POS_prophage_sequence.txt",
                           header = FALSE, sep='#')
colnames(prot_pos_pro) <- c("protein", "coord_start", "coord_stop", "coord_strand", "annotation_info")
prot_pos_pro$protein <- gsub('>', '', prot_pos_pro$protein)
prot_pos_pro$protein <- gsub(' ', '', prot_pos_pro$protein)

pro_gene_annotation <- read_tsv("03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/VIBRANT_annotations_prophage_sequence.tsv")

pro_gene_annotation <- merge(prot_pos_pro, pro_gene_annotation)
colnames(pro_gene_annotation)[c(2,3,4,6)] <- c( 'start', 'end', 'strand', 'seq_id')

pro_virus <- read_paf('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/L34922_LS1/prophage_vs_virus.paf')
pro_virus$seq_id <- 'L34922_LS1'


##### PROPHAGE VS ACTIVE PHAGE

# sequence
my_s0 <- tibble(
  seq_id = c("L34922_LS1", 
             "MGYG000132487_6_RagTag:54887-90884"),
  length = c(34922,35997)
)

# genes
my_g0 <- rbind(L34922_LS1_gene_annotation[,c("seq_id", "start", "end", "strand", "Pfam name", "VOG","VOG name")],
               pro_gene_annotation[,c("seq_id", "start", "end", "strand", "Pfam name", "VOG","VOG name")])
my_g0$strand <- ifelse(my_g0$strand>0, '+', '-')
colnames(my_g0)[ncol(my_g0)] <- 'VOG_name'
my_g0$VOG_name_ed <- sub("^[^ ]* ", "", my_g0$VOG_name)
my_g0$VOG_name_ed <- gsub('hypothetical protein', 'HP', my_g0$VOG_name_ed)
my_g0[!is.na(my_g0$`Pfam name`),]$VOG_name_ed <- my_g0[!is.na(my_g0$`Pfam name`),]$`Pfam name`
my_g0$Color <- "white"
my_g0[!is.na(my_g0$`Pfam name`),]$Color <- 'red'
my_g0[!is.na(my_g0$VOG_name) & is.na(my_g0$`Pfam name`),]$Color <- "blue"
my_g0[my_g0$seq_id=="MGYG000132487_6_RagTag:54887-90884",]$VOG_name_ed <- NA
my_g0[my_g0$seq_id=="LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224",]$seq_id <- "L34922_LS1"


gggenomes(genes=my_g0, seqs=my_s0, links=pro_virus) + 
  geom_seq() +         # draw contig/chromosome lines
  geom_bin_label(size=3) +   # label each sequence 
  geom_gene(aes(fill=Color)) +         # draw genes as arrow
  geom_gene_label(aes(label=VOG_name_ed), size=3, vjust = 1) + 
  geom_link()



##### TWO BIFIDOPHAGES ####

# sequence
my_s1 <- tibble(
  seq_id = c("L34922_LS1", 
             "L37775_LS1"),
  length = c(34922,37775)
)

# genes
my_g1 <- rbind(L34922_LS1_gene_annotation[,c("seq_id", "start", "end", "strand", "Pfam name", "VOG","VOG name")],
               L37775_LS1_gene_annotation[,c("seq_id", "start", "end", "strand", "Pfam name", "VOG","VOG name")])
my_g1$strand <- ifelse(my_g0$strand>0, '+', '-')
colnames(my_g1)[ncol(my_g1)] <- 'VOG_name'
my_g1$VOG_name_ed <- sub("^[^ ]* ", "", my_g1$VOG_name)
my_g1$VOG_name_ed <- gsub('hypothetical protein', 'HP', my_g1$VOG_name_ed)
my_g1[!is.na(my_g1$`Pfam name`),]$VOG_name_ed <- my_g1[!is.na(my_g1$`Pfam name`),]$`Pfam name`
my_g1$Color <- 'Unknown'
my_g1[!is.na(my_g1$`Pfam name`),]$Color <- 'Pfam'
my_g1[!is.na(my_g1$VOG_name) & is.na(my_g1$`Pfam name`),]$Color <- "VOG"
my_g1[my_g1$seq_id=="LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224",]$seq_id <- "L34922_LS1"
my_g1[my_g1$seq_id=="LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438",]$seq_id <- "L37775_LS1"

gggenomes(genes=my_g1[my_g1$seq_id=='L34922_LS1',], seqs=my_s1[my_s1$seq_id=='L34922_LS1',]) + 
  geom_seq() +         # draw contig/chromosome lines
  geom_bin_label(size=3) +   # label each sequence 
  geom_gene(aes(fill=Color)) +         # draw genes as arrow
  geom_gene_label(aes(label=VOG_name_ed), size=3, vjust = 1) +
  scale_fill_manual(values=c('red', 'blue', 'grey')) + 
  theme(legend.position = "bottom")



##### FOR VISUALIZATION ####
write.table(my_g1, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6D_selected_bidiophages_genome_annotation.txt', sep='\t', row.names=F)
write.table(my_s1, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6D_selected_bidiophages_seq_info.txt', sep='\t', row.names=F)
 