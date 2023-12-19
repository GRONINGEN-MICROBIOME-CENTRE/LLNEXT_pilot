setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Supplementary figures
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(ggplot2)
library(ggforce)
library(MetBrewer)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(gggenomes)
library(ggtext)
##############################
# Input data
##############################
###### For Supplementary Figure 1A #####
N_timepoints_MGS <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1A_N_timepoints_MGS.txt", sep = '\t', header=T)

###### For Supplementary Figure 1B #####
N_timepoints_VLP <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1A_N_timepoints_VLP.txt", sep = '\t', header=T)

###### For Supplementary Figure 1C #####
place_delivery_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1B_place_delivery_stat.txt", sep = '\t', header=T)

###### For Supplementary Figure 1D #####
ever_never_BF <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1B_ever_never_BF_stat.txt", sep = '\t', header=T)

###### For Supplementary Figure 1E #####
feeding_p_timepoint <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1B_feeding_mode_p_timepoint_stat.txt", sep = '\t', header=T)
feeding_p_timepoint$Timepoint <- factor(feeding_p_timepoint$Timepoint , levels=c('M1', 'M2', 'M3', 'M6','M9', 'M12'), ordered = T)


###### For Supplementary Figure 1F #####
VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
VLP_metadata$Type <- as.factor(VLP_metadata$Type)

###### For Supplementary Figure 1G #####
bacstability_M1_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1B_N_retained_BacSp_from_M1_infants_stat.txt", sep = '\t', header = T)
bacstability_M1_stat$Timepoint <- factor(bacstability_M1_stat$Timepoint, levels=c('M1', 'M2', 'M3', 'M6','M9', 'M12'), ordered = T)
Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')

bacstability_M1 <- read.table('03a.RESULTS/N_retained_BacSp_from_M1_infants.txt', sep='\t', header=T)
bacstability_M1_dots <- reshape2::melt(bacstability_M1)
bacstability_M1_dots <- bacstability_M1_dots[grep('Not_retained', bacstability_M1_dots$variable, invert=T),]
bacstability_M1_dots$Timepoint <- gsub('.*_','',bacstability_M1_dots$variable)
bacstability_M1_dots$Condition <- 'Retained'
bacstability_M1_dots[grep('Richness', bacstability_M1_dots$variable),]$Condition <- 'Richness'


###### For Supplementary Figure 1H #####
bacstability_M1_abundance_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1C_RA_retained_BacSp_from_M1_infants_stat.txt", sep = '\t', header = T)
bacstability_M1_abundance_stat$Timepoint <- factor(bacstability_M1_abundance_stat$Timepoint, levels=c('M1', 'M2', 'M3', 'M6','M9', 'M12'), ordered = T)

bacstability_M1_abundance <- read.table('03a.RESULTS/Abundance_retained_BacSp_from_M1_infants.txt', sep='\t', header=T)
bacstability_M1_abundance_dots <- reshape2::melt(bacstability_M1_abundance)
bacstability_M1_abundance_dots <- bacstability_M1_abundance_dots[grep('Not_retained', bacstability_M1_abundance_dots$variable, invert = T),]
bacstability_M1_abundance_dots$Timepoint <- gsub('.*_','',bacstability_M1_abundance_dots$variable)
bacstability_M1_abundance_dots$Condition <- 'Retained'
bacstability_M1_abundance_dots[grep('Total_space',bacstability_M1_abundance_dots$variable),]$Condition <- 'Total_space'



###### For Supplementary Figure 1I #####
bacstability_P3_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1D_N_retained_BacSp_from_P3_mothers_stat.txt", sep = '\t', header=T)
bacstability_P3_stat$Timepoint <- factor(bacstability_P3_stat$Timepoint, levels=c("P3","P7", "B",
                                                                                  "M1", "M2", "M3"), ordered = T)
bacstability_P3 <- read.table('03a.RESULTS/N_retained_BacSp_from_P3_mothers.txt', sep='\t', header=T)
bacstability_P3_dots <- reshape2::melt(bacstability_P3)
bacstability_P3_dots <- bacstability_P3_dots[grep('Not_retained', bacstability_P3_dots$variable, invert=T),]
bacstability_P3_dots$Timepoint <- gsub('.*_','',bacstability_P3_dots$variable)
bacstability_P3_dots$Condition <- 'Retained'
bacstability_P3_dots[grep('Richness',bacstability_P3_dots$variable),]$Condition <- 'Richness'


###### For Supplementary Figure 1J #####
bacstability_P3_abundance_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1E_RA_retained_BacSp_from_P3_mothers_stat.txt", sep = '\t', header = T)
bacstability_P3_abundance_stat$Timepoint <- factor(bacstability_P3_abundance_stat$Timepoint, levels=c("P3","P7", "B",
                                                                                                      "M1", "M2", "M3"), ordered = T)
bacstability_P3_abundance <- read.table('03a.RESULTS/Abundance_retained_BacSp_from_P3_mothers.txt', sep='\t', header=T)
bacstability_P3_abundance_dots <- reshape2::melt(bacstability_P3_abundance)
bacstability_P3_abundance_dots <- bacstability_P3_abundance_dots[grep('Not_retained', bacstability_P3_abundance_dots$variable, invert = T),]
bacstability_P3_abundance_dots$Timepoint <- gsub('.*_','',bacstability_P3_abundance_dots$variable)
bacstability_P3_abundance_dots$Condition <- 'Retained'
bacstability_P3_abundance_dots[grep('Total_space',bacstability_P3_abundance_dots$variable),]$Condition <- 'Total_space'


###### For Supplementary Figure 2A #####
Veronese <- met.brewer('Veronese')
p_bac_frac_infants_individual_df <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1FG_microbiome_sorting_df_infants.txt', sep='\t', header=T)
plot_bac_frac_infants_per_sample_N <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1F_N_fractions_personal_microbiome_infants.txt', sep='\t', header=T)
plot_bac_frac_infants_per_sample_N$Infant <- factor(plot_bac_frac_infants_per_sample_N$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered=T)
plot_bac_frac_infants_per_sample_N$Timepoint <- factor(plot_bac_frac_infants_per_sample_N$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M9', 'M12'), ordered=T)

###### For Supplementary Figure 2B #####
p_bac_frac_infants_per_sample_ab <-  read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig1G_RA_fractions_personal_microbiome_infants.txt', sep='\t', header=T)
p_bac_frac_infants_per_sample_ab$Infant <- factor(p_bac_frac_infants_per_sample_ab$Infant, levels=p_bac_frac_infants_individual_df$Infant, ordered = T)
p_bac_frac_infants_per_sample_ab$Timepoint <- factor(p_bac_frac_infants_per_sample_ab$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M9', 'M12'), ordered=T)

###### For Supplementary Figure 2C #####
p_vir_frac_mothers_individual_df <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig2AB_sorting_df_mothers.txt', sep='\t', header=T)
plot_vir_frac_mothers_per_sample_N <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig2A_N_fractions_personal_virome_mothers.txt', sep='\t', header=T)
plot_vir_frac_mothers_per_sample_N$Mother <- factor(plot_vir_frac_mothers_per_sample_N$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered = T)
plot_vir_frac_mothers_per_sample_N$Timepoint <- factor(plot_vir_frac_mothers_per_sample_N$Timepoint, levels=c('P7', 'B', 'M1', 'M2', 'M3'), ordered=T)

###### For Supplementary Figure 2D #####
p_vir_frac_mothers_per_sample_ab <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig2B_RA_fractions_personal_virome_mothers.txt', sep='\t', header=T)
p_vir_frac_mothers_per_sample_ab$Mother <- factor(p_vir_frac_mothers_per_sample_ab$Mother, levels=p_vir_frac_mothers_individual_df$Mother, ordered = T)
p_vir_frac_mothers_per_sample_ab$Timepoint <- factor(p_vir_frac_mothers_per_sample_ab$Timepoint, levels=c('P7', 'B', 'M1', 'M2', 'M3'), ordered=T)

###### For Supplementary Figure 2E #####
plot_bac_frac_mothers_per_sample_N <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig2E_N_fractions_personal_microbiome_mothers.txt', sep='\t', header=T)
p_bac_frac_mothers_individual_df <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig2EF_microbiome_sorting_df_mothers.txt', sep='\t', header=T)
plot_bac_frac_mothers_per_sample_N$Mother <- factor(plot_bac_frac_mothers_per_sample_N$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered = T)
plot_bac_frac_mothers_per_sample_N$Timepoint <- factor(plot_bac_frac_mothers_per_sample_N$Timepoint, levels=c('P3','P7', 'B', 'M1', 'M2', 'M3'), ordered = T)

###### For Supplementary Figure 2F #####
p_bac_frac_mothers_per_sample_ab <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig2F_RA_fractions_personal_microbiome_mothers.txt', sep='\t', header=T)
p_bac_frac_mothers_per_sample_ab$Mother <- factor(p_bac_frac_mothers_per_sample_ab$Mother, levels=p_bac_frac_mothers_individual_df$Mother, ordered = T)
p_bac_frac_mothers_per_sample_ab$Timepoint <- factor(p_bac_frac_mothers_per_sample_ab$Timepoint, levels=c('P3','P7', 'B', 'M1', 'M2', 'M3'), ordered = T)

###### For Supplementary Figure 3A #####
virus_host_interaction_0.2perc <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3A_mixed_model_outcome_vOTUs-host-genus_aggregates_at_0.2_prevalence_vs_bacgen.txt', sep='\t', header = T)
virus_host_interaction_0.2perc$Feature <- paste0("*",virus_host_interaction_0.2perc$Feature,"*")
virus_host_interaction_0.2perc$Feature <- factor(virus_host_interaction_0.2perc$Feature, levels=unique(virus_host_interaction_0.2perc$Feature), ordered = T)
virus_host_interaction_0.2perc$Bug <- gsub('vOTUs_', '', virus_host_interaction_0.2perc$Bug)
virus_host_interaction_0.2perc$Bug <- paste0("*",virus_host_interaction_0.2perc$Bug,"*", ' &phi;')
virus_host_interaction_0.2perc$Bug <- factor(virus_host_interaction_0.2perc$Bug , levels=paste0(unique(virus_host_interaction_0.2perc$Feature), ' &phi;'), ordered = T)

###### For Supplementary Figure 3E #####
PPV_hosts <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Sup_xx_N_overlapping_PP_species_PPV_vs_PPB.txt", sep='\t', header=T)
PPV_hosts$Var1 <- paste0("*", PPV_hosts$Var1, "*")
PPV_hosts[PPV_hosts$Var1=="*Other*",]$Var1 <- "Other"
PPV_hosts$variable <- factor(PPV_hosts$variable, levels=c("N_infecting_vOTUs", "N_overlap_PP_species"), ordered = T)
PPV_hosts$Var1 <- factor(PPV_hosts$Var1, levels=PPV_hosts[1:7,]$Var1, ordered = T)

###### For Supplementary Figure 3B #####
correlations_to_plot <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3B_FDR_correlated_bacteria_vs_its_predicted_phages.txt', sep='\t', header = T)
correlations_to_plot$Timepoint <- factor(correlations_to_plot$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M12'), ordered = T)


###### For Supplementary Figure 3C #####
colorsSFig3AB <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3AB_color_scheme.txt", sep='\t', header=T)
colorsSFig3AB$Genus <- paste0("*", colorsSFig3AB$Genus, "*")

dynamic_vOTU_aggr <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3A_mean_vOTUs_aggregate_genus_p_time_infant.txt", sep='\t', header=T)
dynamic_vOTU_aggr$Timepoint <- factor(dynamic_vOTU_aggr$Timepoint, levels=c('M1', 'M2', 'M3', 'M6','M12'), ordered = T)
dynamic_vOTU_aggr$variable <- paste0("*", dynamic_vOTU_aggr$variable, "*")

###### For Supplementary Figure 3D #####
dynamic_bacgen <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3B_mean_bacterial_genus_p_time_infant.txt", sep='\t', header=T)
dynamic_bacgen$Timepoint <- factor(dynamic_bacgen$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M9', 'M12'), ordered = T)
dynamic_bacgen$variable <- paste0("*", dynamic_bacgen$variable, "*")

###### For Supplementary Figure 4AB #####
MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T)
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("P3","P7", "B","M1", "M2", "M3", 'M6', 'M9', 'M12'), ordered = T)

###### For Supplementary Figure 4C #####
place_of_delivery_bacteria <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3G_A_muciniphila_place_of_birth.txt', sep='\t', header=T)

###### For Supplementary Figure 4D #####
place_of_delivery_phages <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3H_A_muciniphila_phages_place_of_birth.txt', sep='\t', header=T)

###### For Supplementary Figure 4E #####
A_muciniphila_VH <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3I_A_muciniphila_correlation_with_its_phages.txt', sep='\t', header=T)

###### For Supplementary Figure 4F #####
comparison_VLP_VLP_pro <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig4F_increase_in_vOTUs_sharedness_with_prophages.txt", sep='\t', header=T)

###### For Supplementary Figure 4G #####
virus_hists_data <- readRDS("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig4_within_individual_virus_variation_cut_off.rds")
selected_viruses <- read.table('02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', header=T)

###### For Supplementary Figure 5A #####
bacterium_hists_data <- readRDS("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig5_within_individual_bacterium_variation_cut_off.rds")
species_names <- read.table('02.CLEAN_DATA/List_bacterial_strains_reconstructed_and_tested.txt', sep='\t', header=T)
species_names$easy_names <- paste0("*", gsub('_', ' ', species_names$species), "*", "<br>", "(", species_names$Host_SGB, ")")

###### For Supplementary Figure 6A #####
L37775_LS1_to_B.scardovi_patch_ali <- pafr::read_paf("03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/Phages_vs_bacteria_L37775_LS1.paf")
L37775_LS1_to_B.scardovi_patch_ali <- L37775_LS1_to_B.scardovi_patch_ali[L37775_LS1_to_B.scardovi_patch_ali$qname=='C03X03483E05_NODE_35_length_84704_cov_141.800057',]
L37775_LS1_to_B.scardovi_patch_ali <- L37775_LS1_to_B.scardovi_patch_ali[L37775_LS1_to_B.scardovi_patch_ali$tname=='LN_4F04_VL_266_NODE_141_length_37775_cov_36892.259438',]
L37775_LS1_to_B.scardovi_patch_ali$qname <- "B.scardovii patched\ngenome fragment"
L37775_LS1_to_B.scardovi_patch_ali$tname <- "L37775_LS1"

##### for Supplementary Figure 6B #######
L37775_LS1_MGS_reads <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig6B_MGS_read_alignment_patched_B.scardovii.txt", sep='\t', header=T)
L37775_LS1_MGS_reads$color_factor <- factor(L37775_LS1_MGS_reads$color_factor, levels=c("Infant M1", "Infant M2", "Infant M3", "Infant M6", "Infant M9", "Infant M12"), ordered=T)
#L37775__LS1_read_align_colors <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig6", sep='\t', header = T)
#L37775_LS1_MGS_reads$color_factor <- factor(L34922_LS1_MGS_reads$color_factor, levels=c("Other", "Mother P7", "Mother B", 
#                                                                                        "Mother M1", "Mother M2", "Mother M3", 
#                                                                                        "Infant M1", "Infant M2", "Infant M3",
#                                                                                        "Infant M6", "Infant M9"), ordered=T)
L37775_LS1_prophage_region <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig6BC_prophage_region_in_B.scardovii.txt', sep='\t', header=T)
L37775_LS1_MGS_reads <- L37775_LS1_MGS_reads[order(L37775_LS1_MGS_reads$color_factor),]

L37775_LS1_MGS_reads$variable <- factor(L37775_LS1_MGS_reads$variable, levels=unique(L37775_LS1_MGS_reads$variable), ordered=T)

##### for Figure 6F #######
L37775_LS1_VLP_reads <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig6C_VLP_read_alignment_patched_B.scardovii.txt", sep='\t', header=T)
L37775_LS1_VLP_reads$color_factor <- factor(L37775_LS1_VLP_reads$color_factor, levels=c("Infant_M1", "Infant_M2", "Infant_M3", "Infant_M6", "Infant_M12"), 
                                            ordered=T)
L37775_LS1_VLP_reads$variable <- factor(L37775_LS1_VLP_reads$variable, levels=unique(L37775_LS1_VLP_reads$variable), ordered=T)

### Sup 6d
bifidophages_genes <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6D_selected_bidiophages_genome_annotation.txt', sep='\t', header=T)
bifidophages_genes$Color <- factor(bifidophages_genes$Color, levels=c("Pfam", 'VOG', 'Unknown'), ordered = T)
# remove HP:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="HP",]$VOG_name_ed <- NA
# remove predicted but unnanotated:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="Protein D14",]$VOG_name_ed <- NA
# remove predicted but unnanotated:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="Gene 22 protein",]$VOG_name_ed <- NA
# remove predicted but unnanotated:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="Gene 19 protein",]$VOG_name_ed <- NA
# remove predicted but unnanotated:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="Gene 15 protein",]$VOG_name_ed <- NA
# remove predicted but unnanotated:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="Gene 5 protein",]$VOG_name_ed <- NA
# remove predicted but unnanotated:
bifidophages_genes[!is.na(bifidophages_genes$VOG_name_ed) & bifidophages_genes$VOG_name_ed=="gp14",]$VOG_name_ed <- NA

bifidophages_seqs <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6D_selected_bidiophages_seq_info.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################

###### Supplementary Figure 1A #####
SFig1A <- ggplot(data=N_timepoints_MGS, aes(x=N_timepoints, y=N_individuals, fill=Type)) +
  ylim(0,21) +
  geom_bar(stat='identity', position='dodge', alpha=0.5, color="black") +
  geom_vline(aes(xintercept = 5, color="Mode"), linewidth=0.5, lty="dashed") +
  geom_vline(xintercept = 5, color="#00BFC4", linewidth=0.5, lty="dashed") +
  geom_vline(xintercept = 6, color="#F8766D", linewidth=0.5, lty="dashed") +
  labs(x="N timepoints", y="N individuals", fill="Family member", tag="a") +
  scale_color_manual(name = "", values = c(Mode = "grey")) + 
  ggtitle("Total metagenomes") +
  theme_bw() + 
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position="bottom",
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5),
         color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) 

###### Supplementary Figure 1B #####
SFig1B <-  ggplot(data=N_timepoints_VLP, aes(x=N_timepoints, y=N_individuals, fill=Type)) +
  ylim(0,15) +
  geom_bar(stat='identity', position='dodge', alpha=0.5, color="black") +
  geom_vline(aes(xintercept = 4, color="Mode"), linewidth=0.5, lty="dashed") +
  geom_vline(xintercept = 4, color="#00BFC4", linewidth=0.5, lty="dashed") +
  geom_vline(xintercept = 2, color="#F8766D", linewidth=0.5, lty="dashed") +
  labs(x="N timepoints", y="N individuals", fill="Family member", tag="b") +
  scale_color_manual(name = "", values = c(Mode = "grey")) + 
  ggtitle("VLP metaviromes") +
  theme_bw() + 
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position="bottom",
        legend.key.size=unit(0.7, "line")) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5),
         color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) 

###### Supplementary Figure 1C #####
SFig1C <- ggplot(place_delivery_stat, aes(x = "" , y = percentage, fill = fct_inorder(delivery_place))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  labs(fill="Delivery mode and place", tag="c", x = NULL, y = NULL) +
  ggtitle("Delivery mode and place") + 
  scale_fill_manual(values=c("#EEE2DE", "#EA906C", "#B31312")) +
  geom_label_repel(data = place_delivery_stat,
                   aes(y = pos, label = paste0(percentage, "%")),
                   size = 1.75,  show.legend = FALSE) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5)) +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5, size=7), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position="right",
        legend.key.size=unit(0.7, "line"))
 
SFig1ABC <- (SFig1A | SFig1B | SFig1C) + plot_layout(nrow=1, guides = "collect") & 
  theme(legend.position = "bottom", 
        plot.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5)) 

###### Supplementary Figure 1D #####
SFig1D <- ggplot(ever_never_BF, aes(x="", y=percentage, fill=feeding_mode)) +
  labs(fill="Feeding mode", tag="d", x = NULL, y = NULL) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  geom_label_repel(data = ever_never_BF,
                   aes(y = pos, label = paste0(round(percentage, 1), "%")),
                   size = 1.75,  show.legend = FALSE, color="white") +
  theme_void() +
  ggtitle("Feeding mode over time") + 
  scale_fill_manual(values=c("#EA5455", "#002B5B", "#B31312")) +
  theme(plot.title = element_text(size=7, hjust=0.5), 
        plot.tag = element_text(face="bold", size=7, vjust=-7), 
        legend.position="right", 
        legend.title = element_text(size=7),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5))

###### Supplementary Figure 1E #####
SFig1E <- ggplot(feeding_p_timepoint, aes(x=Timepoint, y=value, fill=feeding_mode)) + 
  geom_bar(position = "fill", stat = "identity") + 
  labs(x="Timepoint", y="% infants", fill="Feeding mode", tag="e") + 
  scale_fill_manual(values=c("#7B113A", "#150E56", "#1597BB", "#B31312"),
                    labels=c('Breastfed', 'Formula-fed', "Mixed fed", "NA")) + 
  scale_y_continuous(labels = scales::percent_format()) +
  ggtitle("Feeding mode per timepoint") + 
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = element_text(size=7,hjust=0.5),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.tag = element_text(face="bold", size=7),
        legend.position = "right",
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5))

###### Supplementary Figure 1F #####
SFig1F <- ggplot(VLP_metadata, aes(Timepoint, viral_richness, fill=Type)) + 
  labs (y="Log<sub>10</sub> richness of vOTUs", x="Timepoint", tag="f") + 
  geom_violin(aes(color=Type),alpha=0.1, show.legend = F) + 
  geom_sina(aes(colour=Type), size=0.6, show.legend = F) +
  ggnewscale::new_scale_color() + 
  geom_boxplot(aes(color=Type),width = 0.3, outlier.shape = NA, alpha=0.3, show.legend = F) + 
  scale_color_manual(values=c("#B31312", "#090580")) + 
  ggtitle("vOTUs richness in VLP metaviromes") +
  facet_grid(~Type, scales="free", space = "free") +
  scale_y_continuous(trans='log10') +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title.y = ggtext::element_markdown(size=7),
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        plot.title = element_text(size=7, hjust=0.5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.position="none")


SFig1DEF <- ( SFig1D | SFig1E | SFig1F) + plot_layout(nrow=1, widths = c(3.5,2.5,4), guides = "collect") & 
  theme(legend.position = "bottom", 
        plot.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5))

###### Supplementary Figure 1G #####
SFig1G <-  ggplot(bacstability_M1_stat[bacstability_M1_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_sina(data=bacstability_M1_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=bacstability_M1_stat[bacstability_M1_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  labs(x="Timepoint", y="N detected bacterial species",  fill='Bacterial species', 
       color='Bacterial species', alpha='', tag="g") + 
  ggtitle("Infant bacteriome") +
  theme_bw() + 
  theme(axis.title=element_text(size=7),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from M1', 'All')) + 
  scale_color_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from M1', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained from M1', 'All')) +
  guides(#color="none", 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

###### Supplementary Figure 1H #####
SFig1H <-  ggplot(bacstability_M1_abundance_stat[bacstability_M1_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position="identity") +
  geom_sina(data=bacstability_M1_abundance_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=bacstability_M1_abundance_stat[bacstability_M1_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  labs(x="Timepoint",y="Relative abundance",fill='Bacterial species', color='Bacterial species',alpha='', tag = "h") + 
  ggtitle("Infant bacteriome") +
  theme_bw() + 
  theme(axis.title=element_text(size=7),
        axis.text=element_text(size=5),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        plot.title = element_text(hjust=0.5, size=7),
        legend.key.size=unit(0.7, "line")) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from M1', 'All')) + 
  scale_color_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from M1', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained from M1', 'All')) +
  guides(#color="none",
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top',title.hjust = 0.5), 
         alpha="none")

###### Supplementary Figure 1I #####
SFig1I <- ggplot(bacstability_P3_stat[bacstability_P3_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_sina(data=bacstability_P3_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=bacstability_P3_stat[bacstability_P3_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  labs(x="Timepoint", y="N detected bacterial species",  fill='Bacterial species', 
       color='Bacterial species', alpha='', tag="i") + 
  ggtitle("Mother bacteriome") +
  theme_bw() + 
  theme(axis.title=element_text(size=7),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P3', 'All')) + 
  scale_color_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from P3', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained from P3', 'All')) +
  guides(#color="none", 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

###### Supplementary Figure 1J #####
SFig1J <- ggplot(bacstability_P3_abundance_stat[bacstability_P3_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position="identity") +
  geom_sina(data=bacstability_P3_abundance_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=bacstability_P3_abundance_stat[bacstability_P3_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=1)  +
  labs(x="Timepoint",y="Relative abundance",fill='Bacterial species', color='Bacterial species',alpha='', tag = "j") + 
  ggtitle("Mother bacteriome") +
  theme_bw() + 
  theme(axis.title=element_text(size=7),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P3', 'All')) + 
  scale_color_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from P3', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained from P3', 'All')) +
  guides(#color="none",
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top',title.hjust = 0.5), 
         alpha="none")

SFig1GHIJ <- (SFig1G | SFig1H | SFig1I | SFig1J) + plot_layout(ncol=4, guides = "collect") & 
  theme(legend.position = "bottom", 
        axis.title = element_text(size=7),
        axis.text = element_text(size=5), 
        plot.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5))

SFig1 <- SFig1ABC / SFig1DEF / SFig1GHIJ & theme(legend.margin=margin(0,0,0,0),
                                                 legend.box.margin=margin(-5,-5,-5,-5))

pdf("Figures/Supplementary_figure1/Supplementary_figure_1_1812.pdf", width=18/2.54, height=21/2.54)
SFig1
dev.off()


###### Supplementary Figure 2A #####
SFig2A <-  ggplot(plot_bac_frac_infants_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Infant samples",y="% of sample composition",fill='Fraction', tag="a") +
  ggtitle("Infant bacteriome") +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPB', 'TDB')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'left', title.hjust = 0.5))

###### Supplementary Figure 2B #####
SFig2B <- ggplot(p_bac_frac_infants_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Infant samples",y="% relative abundance",fill='Fraction', tag="b") +
  ggtitle("Infant bacteriome") +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPB', 'TDB')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'left', title.hjust = 0.5))

###### Supplementary Figure 2C #####
SFig2C <- ggplot(plot_vir_frac_mothers_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Maternal samples",y="% of sample composition",fill='Fraction', tag="c") +
  ggtitle("Mother virome") +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPV', 'TDV')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'left', title.hjust = 0.5))

###### Supplementary Figure 2D #####
SFig2D <- ggplot(p_vir_frac_mothers_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Maternal samples",y="% relative abundance",fill='Fraction', tag="d") +
  ggtitle("Mother virome") +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPV', 'TDV')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'left', title.hjust = 0.5))
###### Supplementary Figure 2E #####
SFig2E <- ggplot(plot_bac_frac_mothers_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Maternal samples",y="% of sample composition",fill='Fraction', tag="e") +
  ggtitle("Mother bacteriome") +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPV', 'TDV')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'left', title.hjust = 0.5))

###### Supplementary Figure 2F #####
SFig2F <- ggplot(p_bac_frac_mothers_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Maternal samples",y="% relative abundance",fill='Fraction', tag="f") +
  ggtitle("Mother bacteriome") +
  facet_grid(~Mother, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPV', 'TDV')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_blank(),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'left', title.hjust = 0.5))

SFig2 <- (SFig2A / SFig2B / SFig2C / SFig2D / SFig2E / SFig2F) + plot_layout(guides = "collect") & 
                                                                       theme(legend.position = "top", 
                                                                       legend.text = element_text(size=5),
                                                                       legend.margin=margin(0,0,0,0),
                                                                       legend.box.margin=margin(-5,-5,-5,-5))

pdf("Figures/Supplementary_figure2/Supplementary_figure_2_1812.pdf", width=18/2.54, height=21/2.54)
SFig2
dev.off()

###### Supplementary Figure 3A #####
SFig3A <- ggplot(virus_host_interaction_0.2perc, aes(Feature, Bug, fill = Estimate))+
  labs(x='Bacterial genera', y='vOTUs-host-genus aggregates', tag='a') +
  geom_tile(aes(color=diag), width = 0.98, height = 0.98) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.7,0.7), space = "Lab", 
                       name="Effect size") +
  theme_minimal() + 
  theme(axis.text.x = ggtext::element_markdown(angle = 45, vjust = 1, 
                                   size = 5, hjust = 1),
        axis.text.y = ggtext::element_markdown( vjust = 1, 
                                    size = 5, hjust = 1),
        axis.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        plot.tag = element_text(face="bold", size=7),
        legend.position = "right",
        legend.key.size=unit(0.7, "line")) +
  geom_text(aes(label=Significane_level), color="black", size=1.75) + 
  scale_color_manual(guide = "none", values = c("TRUE" = "black", "FALSE" = "white")) +
  coord_fixed()

###### Supplementary Figure 3B #####
SFig3B <- ggplot(correlations_to_plot, aes(Klebsiella, vOTUs_Klebsiella, color=Timepoint)) +
  geom_point(size=0.7) + 
  labs(x="clr-t abundance of bacterial genus", y="clr-t abundance of vOTUs-aggregate", tag='b') +
  ggtitle("*Klebsiella* and their predicted<br>phages over time") +
  annotate(geom = "text", x = 2.5, y=6, label="FDR=2.2e-06, beta=0.3", size=1.75) +
  geom_smooth(inherit.aes = F, aes(Klebsiella, vOTUs_Klebsiella), method = "lm") +
  theme_bw() + 
  theme(plot.tag = element_text(face="bold", size=7, vjust=-4), 
        plot.title = ggtext::element_markdown(size=7, hjust=0.5),
        axis.text = element_text(size=5),
        axis.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line")) + 
  guides(color=guide_legend(nrow=2, byrow=TRUE, title.position = 'top', title.hjust = 0.5))


###### Supplementary Figure 3C #####
SFig3C <- ggplot(dynamic_vOTU_aggr, aes(x = Timepoint, y = mean, color = variable, group = variable)) +
  geom_point(size=0.7, show.legend = F) +
  geom_line(size=0.5, alpha=0.8, show.legend = F) +
  scale_color_manual(values = colorsSFig3AB$Color) +
  ggtitle("vOTUs-host-genus\naggregates")+
  theme_bw()+
  labs(x="Timepoint", y = "Mean clr-t abundance", color='Bacterial host genera', tag="c")+
  theme(plot.tag = element_text(face="bold", size=7, vjust=-4), 
        plot.title = element_text(size=7, hjust=0.5),
        axis.text = element_text(size=5),
        axis.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = ggtext::element_markdown(size=5),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line")) +
  guides(color=guide_legend(nrow=4, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

###### Supplementary Figure 3D #####
SFig3D <- ggplot(dynamic_bacgen, aes(x = Timepoint, y = mean, color = variable, group = variable)) +
  geom_point(size=0.7) +
  geom_line(size=0.5, alpha=0.8) +
  scale_color_manual(values = colorsSFig3AB$Color) +
  ggtitle("Bacterial genera")+
  theme_bw()+
  labs(x="Timepoint", y = "Mean clr-t abundance", color='Bacterial host genera', tag="d")+
  theme(plot.tag = element_text(face="bold", size=7, vjust=-4), 
        plot.title = ggtext::element_markdown(size=7, hjust=0.5),
        axis.text = element_text(size=5),
        axis.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = ggtext::element_markdown(size=5),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line")) + 
  guides(color=guide_legend(nrow=4, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

###### Supplementary Figure 3E #####
SFig3E <- ggplot(PPV_hosts, aes(x = value, y = Var1, fill = Phylum)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~variable, ncol = 1, scales = 'free_y', drop = TRUE,
             labeller = as_labeller(c(N_overlap_PP_species = "N PPB species overlapping\nwith PPV host prediction", 
                                      N_infecting_vOTUs = "N PPVs predicted to infect\nPPBs at species level")) ) +
  labs(x = "", y = "Bacterial host genera", tag = "e", fill = "Bacterial phylum") +
  #scale_x_log10() +
  theme_bw() +
  theme(axis.text = element_text(size = 5),
        axis.text.x = ggtext::element_markdown(angle=60, vjust=1, hjust = 1, size=5),
        axis.title = element_text(size = 7),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size = 5),
        legend.title = element_text(size=7),
        strip.text = element_text(size=6),
        plot.tag = element_text(face="bold", size=7),
        legend.key.size=unit(0.7, "line")) + 
  coord_flip() + 
  scale_fill_manual(values = c('#FF9C0E', '#1F77B4', '#ED3419', '#B446B3'), name="Bacterial phylum") + 
  guides(fill=guide_legend(nrow=2, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

SFig3BCD <- (SFig3B + SFig3C + SFig3D + SFig3E) + plot_layout(ncol=4, widths = c(2.75, 2,2,2.75), guides="collect") & theme(legend.position = "bottom",
                                                                                                                 legend.margin=margin(0,0,0,0),
                                                                                                                 legend.box.margin=margin(-5,-5,-5,-5))



pdf("Figures/Supplementary_figure3/Supplementary_figure_3_1812.pdf", width=18/2.54, height=21/2.54)
wrap_elements(SFig3A) / wrap_elements(SFig3BCD) + plot_layout(nrow=2, heights = c(6,4)) & theme(legend.position = "bottom",
                                                                                                legend.margin=margin(0,0,0,0),
                                                                                                legend.box.margin=margin(-5,-5,-5,-5))

dev.off()

###### Supplementary Figure 4A #####
SFig4A <- ggplot(MGS_metadata, aes(Timepoint, temperate_RA_MGS, fill=Type)) +
  labs (y="Relative abundance\nof temperate phages", x="Timepoint", tag="a") + 
  geom_violin(aes(color=Type),alpha=0.1, show.legend = F) + 
  geom_sina(aes(colour=Type), size=0.6, show.legend = F) +
  ggnewscale::new_scale_color() + 
  geom_boxplot(aes(color=Type),width = 0.4, outlier.shape = NA, alpha=0.3, show.legend = F, lwd=0.3) + 
  scale_color_manual(values=c("#B31312", "#090580")) + 
  ggtitle("Temperate phages in MGS metaviromes") +
  facet_grid(~Type, scales="free", space = "free") +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        plot.title = element_text(size=7, hjust=0.5), 
        plot.tag = element_text(face="bold", size=7, vjust = -4),
        legend.position="none")

###### Supplementary Figure 4B #####
SFig4B <- ggplot(MGS_metadata[!is.na(MGS_metadata$infant_ever_never_breastfed),], aes(Timepoint, temperate_richness_MGS, fill=infant_ever_never_breastfed)) + 
  labs (y="Log<sub>10</sub> richness<br>of temperate phages", x="Timepoint", tag="b") + 
  geom_boxplot(aes(color=infant_ever_never_breastfed), outlier.shape = NA, alpha=0.5, lwd=0.3) + 
  geom_sina(aes(color=infant_ever_never_breastfed), size=0.6) +
  theme_bw()+
  ggtitle("Temperate phages\nin MGS metaviromes") +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        axis.title.y = ggtext::element_markdown(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5, vjust=-2),
        plot.tag = element_text(face="bold", size=7, vjust = -4)) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Breastfed', 'Exclusively\nformula-fed'),
                    values=c("#EA5455", "#002B5B")) +
  scale_color_manual(name = "Feeding mode", 
                     labels=c('Breastfed', 'Exclusively\nformula-fed'),
                     values=c("#EA5455", "#002B5B")) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5), 
         color=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5))


###### Supplementary Figure 4C #####
SFig4C <- ggplot(place_of_delivery_bacteria, aes(infant_place_delivery, Akkermansia_muciniphila, fill=infant_place_delivery)) + 
  labs (y="clr-t abundance", x="Place of birth", tag='c') + 
  ggtitle("*Akkermansia muciniphila*<br>abundance") +
  ylim(-4, 13) +
  geom_sina(aes(color=infant_place_delivery), size=0.6, show.legend = F) +
  scale_color_manual(name = "Place of birth", 
                     labels=c('Home', 'Hospital'),
                     values=c("#9BCDD2", "#FF8551")) +
  ggnewscale::new_scale_color() +
  geom_boxplot(aes(color=infant_place_delivery), width=0.6, outlier.shape = NA, alpha=0.5, show.legend = F, lwd=0.3) + 
  scale_color_manual(name = "Place of birth", 
                     labels=c('Home', 'Hospital'),
                     values=c("#468B97", "#EF6262")) +
  annotate("text", x=1.5, y=12, label="FDR=0.04, beta hospital=-2.5", size=1.75) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = ggtext::element_markdown(size=7, hjust=0.5, vjust=-2),
        plot.tag = element_text(face="bold", size=7, vjust = -4),
        legend.position = "none",
        legend.key.size=unit(0.7, "line")) +
  scale_fill_manual(name = "Place of birth", 
                    labels=c('Home', 'Hospital'),
                    values=c("#9BCDD2", "#FF8551")) 

###### Supplementary Figure 4D #####
SFig4D <- ggplot(place_of_delivery_phages, aes(infant_place_delivery, Akkermansia.muciniphila, fill=infant_place_delivery)) + 
  labs (y="clr-t abundance", x="Place of birth", tag='d') + 
  ggtitle("*Akkermansia muciniphila* phages") +
  ylim(-4, 11) +
  scale_color_manual(name = "Place of birth", 
                     labels=c('Home', 'Hospital'),
                     values=c("#468B97", "#EF6262")) +
  ggnewscale::new_scale_color() +
  geom_sina(aes(color=infant_place_delivery), size=0.6, show.legend = F) +
  geom_boxplot(aes(color=infant_place_delivery), width=0.5, outlier.shape = NA, alpha=0.5, show.legend = F, lwd=0.3) +
  scale_color_manual(name = "Place of birth", 
                     labels=c('Home', 'Hospital'),
                     values=c("#9BCDD2", "#FF8551")) +
  annotate("text", x=1.5, y=10, label="p-value=0.01, beta hospital=-1.5", size=1.75) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = ggtext::element_markdown(size=7, hjust=0.5),
        plot.tag = element_text(face="bold", size=7, vjust = -4),
        legend.position = "none") +
  scale_fill_manual(name = "Place of birth", 
                    labels=c('Home', 'Hospital'),
                    values=c("#9BCDD2", "#FF8551")) 

###### Supplementary Figure 4E #####
SFig4E <- ggplot(A_muciniphila_VH, aes(Akkermansia_muciniphila, `vOTUs_Akkermansia.muciniphila`, color=infant_place_delivery)) + 
  labs (y="clr-t abundance of\nvOTUs aggregate", x="clr-t abundance of bacterial species", tag='e') + 
  ggtitle("*Akkermansia muciniphila* and its phages") +
  geom_point(size=0.5) +
  geom_smooth(inherit.aes = F, aes(Akkermansia_muciniphila, `vOTUs_Akkermansia.muciniphila`), method="lm") +
  annotate("text", x=4, y=7, label="p-value=5.1e-09, beta bacteria=2.6-01", size=1.75)+
  theme_bw() + 
  scale_color_manual(name = "Place of birth", 
                     labels=c('Home', 'Hospital'),
                     values=c("#9BCDD2", "#FF8551")) +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = ggtext::element_markdown(size=7, hjust=0.5),
        plot.tag = element_text(face="bold", size=7, vjust = -3),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7)) 

###### Supplementary Figure 4F #####
SFig4F <- ggplot(comparison_VLP_VLP_pro, aes(source, value, fill=source)) +
  geom_violin(aes(color=source),alpha=0.3, show.legend = F) +
  labs (y="% Infant vOTUs\nshared with mother", x="Type of data", tag='f') + 
  ylim(0, 1) + 
  geom_sina(aes(color=source), size=0.6, show.legend = F) +
  scale_color_manual(values=c("#F29696", "#4C657E")) +
  ggnewscale::new_scale_color() +
  geom_boxplot(aes(color=source), width=0.3, alpha=0.1, show.legend = F) + 
  scale_color_manual(values=c("#882042", "#36413D")) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        plot.title = ggtext::element_markdown(size=7, hjust=0.5),
        plot.tag = element_text(face="bold", size=7, vjust=-3)) + 
  scale_x_discrete(labels=c('VLP metaviromes','Prophage inclusive\nVLP metaviromes')) +
  annotate(geom="text", x=1.5, y=0.95, 
           label="p-value = 0.001, beta=0.05", size=1.75)

###### Supplementary Figure 4G #####
distance_histograms_all <- list()

for ( virusName in names(virus_hists_data) ) {
  
  distance_histograms_all[[virusName]] <- ggplot(virus_hists_data[[virusName]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(virus_hists_data[[virusName]]$Distance) ), alpha=0.7) + 
    geom_density(aes(x=Distance, fill=Variable), alpha=0.2) +
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=Youden[1], color="Youden_index"), linetype="dashed") + 
    geom_vline(aes(xintercept=FDR_value[1], color="FDR_value"), linetype="dashed") +
    annotate(geom = "text", label=paste0("N=",virus_hists_data[[virusName]]$N_within_comparisons[1]), x=Inf, y=Inf, hjust=+1.1,vjust=+2, size=1.75) +
    ggtitle( selected_viruses$ContigID_easy[match(virusName, selected_viruses$Virus)]  ) +
    theme_bw() + 
    theme(title = element_text(size=7), 
          axis.title = element_text(size=7),
          axis.text = element_text(size=5),
          legend.title = element_text(size=7),
          legend.text = element_text(size=5),
          legend.key.size=unit(0.7, "line"),
          plot.title=element_text(size=7)) +
    scale_color_manual(name="Threshold", 
                       labels=c(Youden_index="Youden index", FDR_value="5% FDR"),
                       values=c(Youden_index="black", FDR_value="red"), 
                       guide = guide_legend(order = 2)) + 
    scale_fill_manual(name="Distribution", 
                      labels=c(Between="Unrelated individual comparison", Within="Same individual comparison"),
                      values=c(Between="#F8766D", Within="#00BFC4"), 
                      guide = guide_legend(order = 1))
  
}


combined_plot <- distance_histograms_all[[1]]

for (h in 2:NROW(virus_hists_data)) {
  combined_plot <- combined_plot + distance_histograms_all[[h]]
}

combined_plot <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "") & theme(legend.position = "bottom")


SFig4ABCD <- (SFig4A + SFig4B + SFig4C) + plot_layout(ncol=3, widths=c(5, 2.5, 2.5),guides="collect") & theme(legend.position = "bottom",
                                                                                                                   legend.margin=margin(0,0,0,0),
                                                                                                                   legend.box.margin=margin(-5,-5,-5,-5))

SFig4DEF <- (SFig4D + SFig4E + SFig4F) + plot_layout(ncol=3,guides="keep") & theme(legend.position = "bottom",
                                                                                                     legend.margin=margin(0,0,0,0),
                                                                                                     legend.box.margin=margin(-5,-5,-5,-5))

SFig4 <- SFig4ABCD / SFig4DEF / combined_plot + plot_layout(nrow=3, heights=c(1.5,1.5,7),guides="keep") & theme(legend.position = "bottom",
                                                                                                                                                   legend.margin=margin(0,0,0,0),
                                                                                                                                                   legend.box.margin=margin(-5,-5,-5,-5))

pdf("Figures/Supplementary_figure4/Supplementary_figure_4_1812.pdf", width=20/2.54, height=27/2.54)
SFig4
dev.off()

###### Supplementary Figure 5 #####
distance_histograms_all <- list()

for (bacteriumName in names(bacterium_hists_data)) {
  
  distance_histograms_all[[bacteriumName]] <- ggplot(bacterium_hists_data[[bacteriumName]], aes(x=Distance, fill=Variable)) + 
    geom_histogram(aes(y = (after_stat(count)/sum(after_stat(count)))*100), position = 'identity', bins=length( unique(bacterium_hists_data[[bacteriumName]]$Distance) ), alpha=0.7) + 
    geom_density(aes(x=Distance, fill=Variable), alpha=0.2) +
    labs(x="Normalized Distance", y="proportion (%)") +
    geom_vline(aes(xintercept=Youden[1], color="Youden_index"), linetype="dashed") + 
    geom_vline(aes(xintercept=FDR_value[1], color="FDR_value"), linetype="dashed") +
    annotate(geom = "text", label=paste0("N=",bacterium_hists_data[[bacteriumName]]$N_within_comparisons[1]), x=Inf, y=Inf, hjust=+1.1,vjust=+2, size=1.75) +
    ggtitle( gsub('_', ' ', species_names$easy_names[match(bacteriumName, species_names$Host_SGB)] ) ) +
    theme_bw() + 
    theme(title = element_text(size=7), 
          axis.title = element_text(size=7),
          axis.text = element_text(size=5),
          legend.title = element_text(size=7),
          legend.text = element_text(size=5),
          plot.title = ggtext::element_markdown(size=6),
          legend.key.size=unit(0.7, "line"),) +
    scale_color_manual(name="Threshold", 
                       labels=c(Youden_index="Youden index", FDR_value="5% FDR"),
                       values=c(Youden_index="black", FDR_value="red"), 
                       guide = guide_legend(order = 2)) + 
    scale_fill_manual(name="Distribution", 
                      labels=c(Between="Unrelated individual comparison", Within="Same individual comparison"),
                      values=c(Between="#F8766D", Within="#00BFC4"), 
                      guide = guide_legend(order = 1))
}


combined_plot <- distance_histograms_all[[1]]

for (h in 2:NROW(bacterium_hists_data)) {
  combined_plot <- combined_plot + distance_histograms_all[[h]]
}

SFig5 <- combined_plot +
  plot_layout(ncol = 5, guides = "collect") + 
  plot_annotation(title = "") & theme(legend.position = "bottom",
                                      legend.margin=margin(0,0,0,0),
                                      legend.box.margin=margin(-5,-5,-5,-5))

pdf("Figures/Supplementary_figure5/Supplementary_figure_5_1812.pdf", width=18/2.54, height=18/2.54)
SFig5
dev.off()

###### Supplementary Figure 6 #####
SFig6A <- pafr::plot_synteny(L37775_LS1_to_B.scardovi_patch_ali, 
                      q_chrom="B.scardovii patched\ngenome fragment",
                      
                      t_chrom="L37775_LS1", 
                      centre=TRUE,
                      xlab="Genome coordinate, kb") + 
  labs(tag="a") +
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size=7), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        plot.tag = element_text(face="bold", size=7)) +
  annotate(geom="text", x = 70000, y=1.5, label="ANI > 99%\ncoverage=100%\ne-value=0", size=1.75) +
  annotate(geom="text", x = 45000, y=0.70, label="Bifidobacterium scardovii patched genome fragment", size=1.75) + 
  annotate(geom="text", x=45000, y=2.3, label="L37775_LS1", size=1.75)


##### Figure 6B #######
SFig6B <- ggplot(L37775_LS1_MGS_reads, aes(w_center/1000, value, group=variable, color=color_factor)) +
  geom_line() + 
  xlim(0, 85) + 
  labs(x="*Bifidobacterium scardovii* patched genome<br>fragment coordinate, kb", y="Log<sub>10</sub> total metagenome mean read depth", color="Timepoint", tag='b') +
  scale_y_log10() + 
  scale_size_manual(values=c(0.5, 0.1), guide=FALSE) +
  #scale_color_manual(values=L34922_LS1_read_align_colors$Color) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.tag = element_text(face="bold", size=7),
        axis.title.x = ggtext::element_markdown(size=7),
        axis.title.y = ggtext::element_markdown(size=7),
        axis.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5)) + 
  guides(color=guide_legend(nrow=2, title.position="top", title.hjust = 0.5)) + 
  ggnewscale::new_scale_color() +
  geom_vline(data=L37775_LS1_prophage_region, aes(xintercept=coordinate/1000, color="Prophage region"), size=0.75, linetype="dashed") + 
  scale_color_manual(name="",
                     values=c(`Prophage region`="black"), 
                     guide = guide_legend(order = 3)) 

##### Figure 6C #######
Fig6C <- ggplot(L37775_LS1_VLP_reads, aes(w_center/1000, value, group=variable, color=color_factor)) +
  geom_line() + 
  xlim(0, 85) + 
  labs(x="*Bifidobacterium scardovii* patched genome<br>fragment coordinate, kb", y="Log<sub>10</sub> VLP metavirome mean read depth", color="Timepoint", tag='c') +
  scale_y_log10() + 
  scale_size_manual(values=c(0.5, 0.1), guide=FALSE) +
  #scale_color_manual(values=L34922_LS1_read_align_colors$Color) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.tag = element_text(face="bold", size=7),
        axis.title.x = ggtext::element_markdown(size=7),
        axis.title.y = ggtext::element_markdown(size=7),
        axis.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5)) + 
  guides(color=guide_legend(nrow=2, title.position="top", title.hjust = 0.5)) + 
  ggnewscale::new_scale_color() +
  geom_vline(data=L37775_LS1_prophage_region, aes(xintercept=coordinate/1000, color="Prophage region"), size=0.75, linetype="dashed") + 
  scale_color_manual(name="",
                     values=c(`Prophage region`="black"), 
                     guide = guide_legend(order = 3)) 

SFig6D <-gggenomes(genes=bifidophages_genes[bifidophages_genes$seq_id=='L37775_LS1',], 
                    seqs=bifidophages_seqs[bifidophages_seqs$seq_id=='L37775_LS1',]) + 
  geom_seq() +         # draw contig/chromosome lines
  geom_bin_label(size=1.75) +   # label each sequence 
  geom_gene(aes(fill=Color)) +         # draw genes as arrow
  geom_gene_label(aes(label=VOG_name_ed), size=1.75, vjust = 1) +
  labs(tag='d', fill='Genome coordinate, kb\nGene annotation source') +
  scale_fill_manual(values=c('#205375', '#EFEFEF', '#F66B0E')) + 
  theme(legend.position = "bottom",
        plot.tag=element_text(face="bold", size=7, vjust=-3.5), 
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        axis.title = element_text(size=7),
        axis.text = element_text(size=5),
        plot.title = element_text(size=7),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(title.position = "top", title.hjust=0.5))


pdf("Figures/Supplementary_figure6/Supplementary_figure_6_1812.pdf", width=18/2.54, height=17/2.54)
(SFig6A | SFig6B | Fig6C) / SFig6D & theme(legend.position = "bottom",
                                          legend.margin=margin(0,0,0,0),
                                          legend.box.margin=margin(-5,-5,-5,-5))
dev.off()

##############################
# OUTPUT
##############################
