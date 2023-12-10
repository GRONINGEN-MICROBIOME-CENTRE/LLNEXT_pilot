#### LIBRARIES ####
library(tidyverse)

#### INPUT ####
table_of_origin <- read.table('/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination/table_of_origin_noneg405_99_der95', sep = '\t', header=T)
table_of_origin$length = sapply(strsplit(table_of_origin$V1, "_"), function(x) x[8])

rbs_J_50_90 <- read.table('/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination/RbS_proteins_only_50_90', sep='', header = F)
colnames(rbs_J_50_90) <- c("N_rbs_prot_J_50_90", "V1")

vcontact2_output <- read.table("/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination/pilot_crAss_ProkaryoticViralRefSeq211.VC.txt", sep='\t', header=T)
colnames(vcontact2_output)[2] <- "V1"
vcontact2_output <- subset(vcontact2_output, vcontact2_output$V1 %in% table_of_origin$V1)

SILVA_rrna_contigs <- read.table('/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination/PILOT_vir_rrna_contigs.ids')
SILVA_rrna_contigs$vir_rrna_contigs <- T

#### PROCESSING ####
decontamination_df <- 
								list(table_of_origin, vcontact2_output, rbs_J_50_90, SILVA_rrna_contigs) %>% 
								reduce(full_join, by='V1')
decontamination_df[is.na(decontamination_df$N_rbs_prot_J_50_90),]$N_rbs_prot_J_50_90 <- 0
decontamination_df[!is.na(decontamination_df$vir_rrna_contigs),]$vir_rrna_contigs <- 1
decontamination_df[is.na(decontamination_df$vir_rrna_contigs),]$vir_rrna_contigs <- 0

decontamination_df$excluding <- F

#exclude if >=1 ribosomal protein gene and < 3 pVOGs per 10 kb length, VirSorter-negative and non-circular
decontamination_df$excluding[decontamination_df$N_rbs_prot_J_50_90 >=1 & decontamination_df$N_pVOGs_per_10kb < 3 & (decontamination_df$VirSorter==0 & decontamination_df$Topology==0)] <- T
#exclude if having > 3 ribosomal protein genes were excluded from the dataset
decontamination_df$excluding[decontamination_df$N_rbs_prot_J_50_90 > 3] <- T
#exclude if rRNA-positive:
decontamination_df$excluding[decontamination_df$vir_rrna_contigs==T] <- T

vc_to_exclude <- aggregate(decontamination_df$excluding[!is.na(decontamination_df$VC)], by=list(decontamination_df$VC[!is.na(decontamination_df$VC)]), FUN=max)
colnames(vc_to_exclude)[c(1,2)] <- c("VC", "final_exclusion")
decontamination_df <- merge(decontamination_df, vc_to_exclude, by = "VC", all = TRUE)
decontamination_df$final_exclusion[is.na(decontamination_df$final_exclusion)] <- ifelse((is.na(decontamination_df$final_exclusion[is.na(decontamination_df$final_exclusion)]) & decontamination_df$excluding[is.na(decontamination_df$final_exclusion)] == T), 1, 0)
decontamination_df$final_exclusion <- as.logical(decontamination_df$final_exclusion)

#keep if were circular and having >= 1 pVOGs
decontamination_df$final_exclusion[decontamination_df$Topology==T & decontamination_df$N_pVOGs_per_10kb >= 1] <- F
#keep if circular and VirSorter-positive
decontamination_df$final_exclusion[decontamination_df$Topology==T & decontamination_df$VirSorter==T] <- F
#keep if VirSorter-positive and having no ribosomal protein genes
decontamination_df$final_exclusion[decontamination_df$VirSorter==T & decontamination_df$N_rbs_prot_J_50_90 ==0 ] <- F

table_of_origin_clean <- subset(decontamination_df, decontamination_df$final_exclusion==F)

#### OUTPUT ####
write.table(table_of_origin_clean, "/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination/table_of_origin_decontaminated", sep='\t', quote=F, row.names=F)
write.table(table_of_origin_clean$V1, "/scratch/umcg-sgarmaeva/ANALYSIS_NEXT_PILOT/clean_viral/decontamination/viral_IDs_noneg405_99_der95_decontaminated", quote=F, row.names=F, col.names=F)
