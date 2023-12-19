setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

##############################
# Loading libraries
##############################
library(ggplot2)
library(ggrepel)
library(ggExtra)
library(ggforce)
library(patchwork)
library(tidyverse)
library(MetBrewer)
library(EnhancedVolcano)
library(ggforestplot)
library(corrplot)
library(ggtree)
library(pafr)
library(gggenomes)
library(ggtext)
##############################
# Input data
##############################
VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
VLP_metadata$Type <- as.factor(VLP_metadata$Type)
VLP_metadata$Short_sample_ID <- row.names(VLP_metadata)

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T, row.names = "Short_sample_ID")
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("P3", "P7", "B",
                                                                    "M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
MGS_metadata$Type <- as.factor(MGS_metadata$Type)
MGS_metadata$Short_sample_ID <- row.names(MGS_metadata)

##### for Figure 1B #######
vOTUs_data.scores <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1B_NMDS_scores_vOTUs.txt", sep='\t', header=T)
# in order to be able to merge the legends for Fig1B and Fig1C together, we need to simulate one M9 point that will not 
# affect anything except legend
vOTUs_data.scores["C090000V",] <- NA
vOTUs_data.scores["C090000V","Type"] <- "Infant"
vOTUs_data.scores["C090000V","Timepoint"] <- "M9"
# we give it exact coordinates of one of the M12 points so that it will be hidden behind it
vOTUs_data.scores["C090000V","NMDS1"] <- 0.118678861
vOTUs_data.scores["C090000V","NMDS2"] <- 0.2251552
vOTUs_data.scores$Timepoint <- factor(vOTUs_data.scores$Timepoint, levels = c('M1', 'M2', 'M3', 'M6', 'M9', 'M12', 'Mother'), ordered = T)

vOTUs_centroids <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1B_NMDS_centroids_vOTUs.txt", sep='\t', header=T)
vOTUs_centroids$Timepoint <- factor(vOTUs_centroids$Timepoint, levels = c('M1', 'M2', 'M3', 'M6','M12', 'Mother'), ordered = T)
vOTUs_spp.scrs <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1B_NMDS_vectors_vOTUs.txt", sep='\t', header=T)

##### for Figure 1C #######
BacSp_data.scores <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1C_NMDS_scores_BacSp.txt", sep='\t', header=T)
BacSp_data.scores$Timepoint <- factor(BacSp_data.scores$Timepoint, levels = c('M1', 'M2', 'M3', 'M6', 'M9', 'M12', 'Mother'), ordered = T)

BacSp_centroids <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1C_NMDS_centroids_BacSp.txt", sep='\t', header=T)

BacSp_spp.scrs <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1C_NMDS_vectors_BacSp.txt", sep='\t', header=T)

##### for Figure 1D #######
vOTUs_data.scores_mothers <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS_scores_vOTUs_mothers.txt", sep='\t', header=T)
# in order to be able to merge the legends for Fig1D and Fig1E together, we need to simulate one P3 point that will not 
# affect anything except legend
vOTUs_data.scores_mothers["AP30000V",] <- NA
vOTUs_data.scores_mothers["AP30000V","Type"] <- "Mother"
vOTUs_data.scores_mothers["AP30000V","Timepoint"] <- "P3"
# we give it exact coordinates of one of the P7 points so that it will be hidden behind it
vOTUs_data.scores_mothers["AP30000V","NMDS1"] <- -0.3134314794
vOTUs_data.scores_mothers["AP30000V","NMDS2"] <- -0.017969975
vOTUs_data.scores_mothers$Timepoint <- factor(vOTUs_data.scores_mothers$Timepoint, levels = c('P3', 'P7','B','M1', 'M2', 'M3'), ordered = T)

vOTUs_centroids_mothers <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS_centroids_vOTUs_mothers.txt", sep='\t', header=T)
vOTUs_centroids_mothers$Timepoint <- factor(vOTUs_centroids_mothers$Timepoint, levels = c('P7','B','M1', 'M2', 'M3'), ordered = T)
vOTUs_spp.scrs_mothers <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS_vectors_vOTUs_mothers.txt", sep='\t', header=T)

Peru1 <- met.brewer('Peru1')

##### for Figure 1E #######
BacSp_data.scores_mothers <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS_scores_BacSp_mothers.txt", sep='\t', header=T)
BacSp_data.scores_mothers$Timepoint <- factor(BacSp_data.scores_mothers$Timepoint, levels = c('P3', 'P7','B','M1', 'M2', 'M3'), ordered = T)
BacSp_centroids_mothers <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS_centroids_BacSp_mothers.txt", sep='\t', header=T)
BacSp_spp.scrs_mothers <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS_vectors_BacSp_mothers.txt", sep='\t', header=T)

##### for Figure 2A #######
significance_annotation_F2A <- data.frame(label = c("p-value=4.1e-4, beta=4.7e-03", "p-value=1, beta=4e-04"),
                                            Type   = c("Infant", "Mother") )

##### for Figure 2B #######
significance_annotation_F2B <- data.frame(label = c("p-value=1.6e-12, beta=2.7e-03", "p-value=0.3, beta=1.1e-02"),
                                          Type   = c("Infant", "Mother") )

##### for Figure 2C #######
Demuth <- met.brewer("Demuth")
Monet <- met.brewer('Monet')
virstability_M1_stat <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2A_N_retained_vOTUs_from_M1_infants_stat.txt', sep='\t', header=T)
virstability_M1_stat$Timepoint <- factor(virstability_M1_stat$Timepoint, levels = c('M1', 'M2', 'M3', 'M6', 'M12'), ordered = T)

virstability_M1 <- read.table('03a.RESULTS/N_retained_vOTUs_from_M1_infants.txt', sep='\t', header=T)
virstability_M1_dots <- reshape2::melt(virstability_M1)
virstability_M1_dots <- virstability_M1_dots[grep('Not_retained', virstability_M1_dots$variable, invert=T),]
virstability_M1_dots$Timepoint <- gsub('.*_','',virstability_M1_dots$variable)
virstability_M1_dots$Condition <- 'Retained'
virstability_M1_dots[grep('Richness',virstability_M1_dots$variable),]$Condition <- 'Richness'

##### for Figure 2D #######
virstability_M1_abundance_stat <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2B_RA_retained_vOTUs_from_M1_infants_stat.txt', sep='\t', header=T)
virstability_M1_abundance_stat$Timepoint <- factor(virstability_M1_abundance_stat$Timepoint, levels = c('M1', 'M2', 'M3', 'M6', 'M12'), ordered = T)
virstability_M1_abundance <- read.table('03a.RESULTS/Abundance_retained_vOTUs_from_M1_infants.txt', sep='\t', header=T)
virstability_M1_abundance_dots <- reshape2::melt(virstability_M1_abundance)
virstability_M1_abundance_dots <- virstability_M1_abundance_dots[grep('Not_retained', virstability_M1_abundance_dots$variable, invert = T),]
virstability_M1_abundance_dots$Timepoint <- gsub('.*_','',virstability_M1_abundance_dots$variable)
virstability_M1_abundance_dots$Condition <- 'Retained'
virstability_M1_abundance_dots[grep('Total_space',virstability_M1_abundance_dots$variable),]$Condition <- 'Total_space'

##### for Figure 2E #######
virstability_P7_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2C_N_retained_vOTUs_from_P7_mothers_stat.txt", sep='\t', header=T)
virstability_P7_stat$Timepoint <- factor(virstability_P7_stat$Timepoint, levels=c("P7", "B",
                                                                                  "M1", "M2", "M3"), ordered = T)
virstability_P7 <- read.table('03a.RESULTS/N_retained_vOTUs_from_P7_mothers.txt', sep='\t', header=T)
virstability_P7_dots <- reshape2::melt(virstability_P7)
virstability_P7_dots <- virstability_P7_dots[grep('Not_retained', virstability_P7_dots$variable, invert=T),]
virstability_P7_dots$Timepoint <- gsub('.*_','',virstability_P7_dots$variable)
virstability_P7_dots$Condition <- 'Retained'
virstability_P7_dots[grep('Richness',virstability_P7_dots$variable),]$Condition <- 'Richness'

##### for Figure 2F #######
virstability_P7_abundance_stat <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2D_RA_retained_vOTUs_from_P7_mothers_stat.txt', sep='\t', header=T)
virstability_P7_abundance_stat$Timepoint <- factor(virstability_P7_stat$Timepoint, levels=c("P7", "B",
                                                                                            "M1", "M2", "M3"), ordered = T)
virstability_P7_abundance <- read.table('03a.RESULTS/Abundance_retained_vOTUs_from_P7_mothers.txt', sep='\t', header=T)
virstability_P7_abundance_dots <- reshape2::melt(virstability_P7_abundance)
virstability_P7_abundance_dots <- virstability_P7_abundance_dots[grep('Not_retained', virstability_P7_abundance_dots$variable, invert = T),]
virstability_P7_abundance_dots$Timepoint <- gsub('.*_','',virstability_P7_abundance_dots$variable)
virstability_P7_abundance_dots$Condition <- 'Retained'
virstability_P7_abundance_dots[grep('Total_space',virstability_P7_abundance_dots$variable),]$Condition <- 'Total_space'

##### for Figure 2G #######
Veronese <- met.brewer('Veronese')

p_vir_frac_infants_individual_df <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2BD_sorting_df_infants.txt', sep = '\t', header=T)
p_vir_frac_infants_individual_df$facet_name <- paste0('I', sprintf("%02s", 1:14))
plot_vir_frac_infants_per_sample_N <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2C_N_fractions_personal_virome_infants.txt', sep='\t', header=T)
plot_vir_frac_infants_per_sample_N$Infant <- p_vir_frac_infants_individual_df$facet_name[match(plot_vir_frac_infants_per_sample_N$Infant, p_vir_frac_infants_individual_df$Infant)]
plot_vir_frac_infants_per_sample_N$Infant <- factor(plot_vir_frac_infants_per_sample_N$Infant, levels=p_vir_frac_infants_individual_df$facet_name, ordered = T)
plot_vir_frac_infants_per_sample_N$Timepoint <- factor(plot_vir_frac_infants_per_sample_N$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M12'), ordered=T)

##### for Figure 2H #######
p_vir_frac_infants_per_sample_ab <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2D_RA_fractions_personal_virome_infants.txt', sep='\t', header=T)
p_vir_frac_infants_per_sample_ab$Infant <- p_vir_frac_infants_individual_df$facet_name[match(p_vir_frac_infants_per_sample_ab$Infant, p_vir_frac_infants_individual_df$Infant)]
p_vir_frac_infants_per_sample_ab$Infant <- factor(p_vir_frac_infants_per_sample_ab$Infant, levels=p_vir_frac_infants_individual_df$facet_name, ordered = T)
p_vir_frac_infants_per_sample_ab$Timepoint <- factor(p_vir_frac_infants_per_sample_ab$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M12'), ordered=T)

##### for Figure 2I #######
HostGen_average <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2G_HostGen_vOTUs_average_abundance_over_time.txt', sep = '\t', header=T)
HostGen_average$variable <- factor(HostGen_average$variable, levels = c("P3", "P7", "B",
                                                                        "M1", "M2", "M3",
                                                                        "M6", "M9", "M12"), ordered = T)
##### for Figure 2J #######
BacGen_average <-  read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2H_BacGen_average_abundance_over_time.txt', sep='\t', header=T)
BacGen_average$variable <- factor(BacGen_average$variable, levels = c("P3", "P7", "B",
                                                                      "M1", "M2", "M3",
                                                                      "M6", "M9", "M12"), ordered = T)

Gen_palette <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2GH_gen_palette_colors.txt', sep='\t', header=T, as.is = T)

# to further merge the legends: 
simulate_missing_bacgen <- data.frame( setdiff(Gen_palette$Genus, unique(HostGen_average$taxa)), 
                                       rep('P3', length(setdiff(Gen_palette$Genus, unique(HostGen_average$taxa)))),
                                       rep(0, length(setdiff(Gen_palette$Genus, unique(HostGen_average$taxa)))),
                                       rep('Mother', length(setdiff(Gen_palette$Genus, unique(HostGen_average$taxa)))),
                                       rep('MGS', length(setdiff(Gen_palette$Genus, unique(HostGen_average$taxa)))))
colnames(simulate_missing_bacgen) <- colnames(HostGen_average)

HostGen_average <- rbind(HostGen_average, simulate_missing_bacgen)
HostGen_average$Sequencing <- ifelse(HostGen_average$Sequencing=="VLP", "Metavirome", "Whole metavirome")

simulate_missing_hostgen <- data.frame( setdiff(Gen_palette$Genus, unique(BacGen_average$taxa)), 
                                        rep('P3', length(setdiff(Gen_palette$Genus, unique(BacGen_average$taxa)))),
                                        rep(0, length(setdiff(Gen_palette$Genus, unique(BacGen_average$taxa)))),
                                        rep('Mother', length(setdiff(Gen_palette$Genus, unique(BacGen_average$taxa)))))
colnames(simulate_missing_hostgen) <- colnames(BacGen_average)

BacGen_average <- rbind(BacGen_average, simulate_missing_hostgen)

##### for Figure 3A #######
significance_annotation_F3A <- data.frame(label = c("p-value=2.9e-04 beta=-8.7e-02", "p-value=9.2e-02 beta=-5.3e-01"),
                                           Type   = c("Infant", "Mother") )


##### for Figure 3D #######
diff_ab_temperate <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig3C_Differentially_prevalent_aggregated_temperates_Mixed_Models_Feeding.txt", sep='\t', header=T)

##### for Figure 3E #######
sharedness_timepoints_melt <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig3D_infant_vOTUs_sharedness_with_pre-birth_vs_post-birth_samples.txt", sep='\t', header=T)
sharedness_timepoints_melt$Timepoint <- factor(sharedness_timepoints_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M12'), ordered = T)
sharedness_timepoints_melt$variable <- factor(sharedness_timepoints_melt$variable, levels = c('perc_to_pre','perc_to_post'), ordered=T)

##### for Figure 3F #######
sharedness_timepoints_iVM_mVM_melt <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig3E_infant_vOTUs_sharedness_with_pre-birth_vs_post-birth_samples_plus_prophages", sep = '\t', header=T)
sharedness_timepoints_iVM_mVM_melt$Timepoint <- factor(sharedness_timepoints_iVM_mVM_melt$Timepoint, levels=c('M1', 'M2', 'M3', 'M6', 'M12'), ordered = T)
sharedness_timepoints_iVM_mVM_melt$variable <- factor(sharedness_timepoints_iVM_mVM_melt$variable, levels = c('perc_to_pre','perc_to_post'), ordered=T)


##### for Figure 4A #######
plot_distances_select <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/4A_Virus_distances_related_vs_unrelated_FDR_significant.txt', sep='\t', header = T)
N_pairwise_distance <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/4A_Virus_N_pairwise_distance.txt', sep='\t', header = T)

##### for Figure 4B #######
for_plot <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/4B_Perc_transmitted_viruses_in_pairs_maximized_Youden_with_N.txt', sep='\t', header=T)

##### for Figure 4C #######
plot_distances_select_bac <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig4C_Bacterium_distances_related_vs_unrelated.txt', sep='\t', header=T)
plot_distances_select_bac$ord <- sprintf('%02s', plot_distances_select_bac$ord)
N_pairwise_distance_bac <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig4C_N_pairwise_distance_related_vs_unrelated.txt', sep='\t', header=T)
N_pairwise_distance_bac$ord <- sprintf('%02s', N_pairwise_distance_bac$ord)

##### for Figure 4D #######
for_plot_bac <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig4D_Perc_transmitted_bacteria_in_pairs_maximized_Youden_with_N.txt', sep='\t', header = T)
for_plot_bac$ord <- sprintf('%02s', for_plot_bac$ord)

##### for Figure 5A #######
bgcolors <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5A_bg_colors_bacteria_virus_pair.txt", sep='\t', header=T)
bgcolors <- as.matrix(bgcolors)

co_transmission_cor <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5A_cotransmission_correlation_matrix.txt", sep='\t', header=T)
co_transmission_cor <- as.matrix(co_transmission_cor)
row.names(co_transmission_cor) <- gsub("^(.*?)_", "\\1 ", row.names(co_transmission_cor))

co_transmission_FDR <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5A_cotransmission_FDR_matrix.txt", sep='\t', header=T)
co_transmission_FDR <- as.matrix(co_transmission_FDR)
row.names(co_transmission_FDR) <- gsub("^(.*?)_", "\\1 ", row.names(co_transmission_FDR))

##### for Figure 5B #######
cotransmission_stat <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5B_cotransmission_stat.txt", sep='\t', header=T)
Kandinsky <- met.brewer("Kandinsky")

##### for Figure 5C #######
phylo_L_85266_LS0 <- readRDS(file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_L_85266_LS0.rds")
phylo_Host_L_85266_LS0 <- readRDS(file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_Host_L_85266_LS0.rds")
L_85266_LS0_label_colors <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_L_85266_LS0_label_colors.txt", sep='\t', header=T)
L_85266_LS0_label_colors_family <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_L_85266_LS0_family_colors.txt", sep='\t', header=T)
Host_L_85266_LS0_label_colors_bac <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_Host_L_85266_LS0_label_colors.txt", sep='\t', header=T)

##### for Figure 6A #######
phylo_L34922_LS1 <- readRDS(file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_L34922_LS1.rds")
L34922_LS1_label_colors <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_L34922_LS1_label_colors.txt", sep='\t', header=T)
L34922_LS1_label_colors_family <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_L34922_LS1_label_family_colors.txt", sep='\t', header=T)

phylo_Host_L34922_LS1 <- readRDS(file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_Host_L34922_LS1_concurrent.rds")
Host_L34922_LS1_label_colors_bac <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_Host_L34922_LS1_label_colors.txt", sep='\t', header=T)

##### for Figure 6B #######
Host_L34922_LS1_copresence <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6B_L34922_LS1_copresence.txt", sep='\t', header=T)
Host_L34922_LS1_copresence$Timepoint <- factor(Host_L34922_LS1_copresence$Timepoint, levels = c("P3", "P7", "B",
                                                                                                "M1", "M2", "M3",
                                                                                                "M6", "M9", "M12"), ordered = T)
Host_L34922_LS1_copresence$VLP <- as.factor(Host_L34922_LS1_copresence$VLP)
Host_L34922_LS1_copresence$MGS <- as.factor(Host_L34922_LS1_copresence$MGS)
##### for Figure 6C #######
L34922_LS1_to_B.bifidum_patch_ali <- pafr::read_paf("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6C_L34922_LS1_to_B.bifidum_patch.paf")
L34922_LS1_to_B.bifidum_patch_ali$qname <- 'L34922_LS1'
L34922_LS1_to_B.bifidum_patch_ali[L34922_LS1_to_B.bifidum_patch_ali$tname=='MGYG000132487_6_RagTag',]$tname <- 'B. bifidum patched\ngenome fragment'
L34922_LS1_to_B.bifidum_patch_ali <- L34922_LS1_to_B.bifidum_patch_ali[L34922_LS1_to_B.bifidum_patch_ali$nmatch > 55,]

##### for Figure 6D #######
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



##### for Figure 6E #######
L34922_LS1_MGS_reads <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6E_MGS_read_alignment_patched_B.bifidum.txt", sep='\t', header=T)
L34922_LS1_read_align_colors <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6EF_read_alignment_colors.txt", sep='\t', header = T)
L34922_LS1_MGS_reads$color_factor <- factor(L34922_LS1_MGS_reads$color_factor, levels=c("Other", "Mother P7", "Mother B", 
                                                                                        "Mother M1", "Mother M2", "Mother M3", 
                                                                                        "Infant M1", "Infant M2", "Infant M3",
                                                                                        "Infant M6", "Infant M9"), ordered=T)
L34922_LS1_prophage_region <- read.table('02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6EF_prophage_region_in_B.bifidum.txt', sep='\t', header=T)
L34922_LS1_MGS_reads <- L34922_LS1_MGS_reads[order(L34922_LS1_MGS_reads$color_factor),]

L34922_LS1_MGS_reads$variable <- factor(L34922_LS1_MGS_reads$variable, levels=unique(L34922_LS1_MGS_reads$variable), ordered=T)

##### for Figure 6F #######
L34922_LS1_VLP_reads <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6F_VLP_read_alignment_patched_B.bifidum.txt", sep='\t', header=T)
L34922_LS1_VLP_reads$color_factor <- factor(L34922_LS1_VLP_reads$color_factor, levels=c("Other", "Mother_P7", "Mother_B", 
                                                                                        "Mother_M1", "Mother_M2", "Mother_M3", 
                                                                                        "Infant_M2", "Infant_M3"), 
                                            ordered=T)
L34922_LS1_VLP_reads$variable <- factor(L34922_LS1_VLP_reads$variable, levels=unique(L34922_LS1_VLP_reads$variable), ordered=T)


##############################
# ANALYSIS
##############################

###### Figure 1 ########

##### Figure 1B #######
pv <- ggplot(data = vOTUs_data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 1.5, alpha=0.8) + 
  geom_point(data=vOTUs_centroids, aes(fill=Timepoint),shape=23, size=2, color='black') + 
  geom_segment(data = vOTUs_spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = vOTUs_spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  #scale_fill_manual(values = c("#3E0751", "#423A7F", "#3F678B", "#468E8B","#5FB57E", "#9FD55C", "#F9E855") ) +
  xlim(-0.8, 0.8)+
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size = 7), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

xAxisBoxPlot_v <- ggplot(data = vOTUs_data.scores[vOTUs_data.scores$Timepoint!="M9",], aes(x = Timepoint, y = NMDS1)) +
  geom_boxplot(aes(fill=Timepoint), alpha=0.7, show.legend = FALSE) + 
  ylim(-0.8,0.8)+
  labs(tag = 'b') +
  coord_flip() +
  ggtitle("Shift in infant virome composition over time") +
  #scale_fill_manual(values = c("#3E0751", "#423A7F", "#3F678B", "#468E8B","#9FD55C", "#F9E855") ) +
  theme_void() +
  theme(plot.title = element_text(hjust=0.0, size=7), plot.tag = element_text(face="bold", vjust=-7, size=7))

Fig1B <- xAxisBoxPlot_v / pv + 
  plot_layout(heights = c(2,8)) & guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         alpha="none")


##### Figure 1C #######
pb <- ggplot(data = BacSp_data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 1.5, alpha=0.8) + 
  geom_point(data=BacSp_centroids, aes(fill=Timepoint),shape=23, size=2, color='black', ) + 
  geom_segment(data = BacSp_spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = BacSp_spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  xlim(-1.8, 2.55)+
  #scale_fill_manual(values = c("#3E0751", "#423A7F", "#3F678B", "#468E8B", "#5FB57E","#9FD55C", "#F9E855")) +
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size = 7), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position = "bottom") 

xAxisBoxPlot_b <- ggplot(data = BacSp_data.scores, aes(x = Timepoint, y = NMDS1)) +
  geom_boxplot(aes(fill=Timepoint), alpha=0.7, show.legend = FALSE) + 
  ylim(-1.8, 2.55)+
  ggtitle("Shift in infant bacteriome composition over time") +
  labs(tag = 'c') +
  coord_flip() +
  theme_void() +
  theme(plot.title = element_text(hjust=0.0, size=7), plot.tag = element_text(face="bold", vjust=-7, size=7))

Fig1C <- xAxisBoxPlot_b / pb + 
  plot_layout(heights = c(2,8)) + guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         alpha="none")

##### Figure 1D #######
pvm <- ggplot(data = vOTUs_data.scores_mothers, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 1.5, alpha=0.8) + 
  geom_point(data=vOTUs_centroids_mothers, aes(fill=Timepoint),shape=23, size=2, color='black') + 
  geom_segment(data = vOTUs_spp.scrs_mothers,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_label_repel(data = vOTUs_spp.scrs_mothers, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  scale_fill_manual(values=Peru1) +
  scale_color_manual(values=Peru1) +
  xlim(-1,1.05)+
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size = 7), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position = "bottom") +
  guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

xAxisBoxPlot_vm <- ggplot(data = vOTUs_data.scores_mothers[vOTUs_data.scores_mothers$Timepoint!="P3",], aes(x = Timepoint, y = NMDS1)) +
  geom_boxplot(aes(fill=Timepoint), alpha=0.7, show.legend = FALSE) + 
  labs(tag = 'd') +
  ylim(-1,1.05)+
  coord_flip() +
  ggtitle("Shift in mother virome composition over time") +
  scale_fill_manual(values=Peru1[2:6]) +
  theme_void() +
  theme(plot.title = element_text(hjust=0.0, size=7), plot.tag = element_text(face="bold", vjust=-7, size=7))

Fig1D <- xAxisBoxPlot_vm / pvm + 
  plot_layout(heights = c(2,8)) & guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         alpha="none")

##### Figure 1E #######
pbm <- ggplot(data = BacSp_data.scores_mothers, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 1.5, alpha=0.8) + 
  geom_point(data=BacSp_centroids_mothers, aes(fill=Timepoint),shape=23, size=2, color='black', ) + 
  geom_segment(data = BacSp_spp.scrs_mothers,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = BacSp_spp.scrs_mothers, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  scale_fill_manual(values=Peru1) +
  scale_color_manual(values=Peru1) +
  xlim(-1,1)+
  theme_bw()+
  theme(axis.text=element_text(size=5),
        axis.title = element_text(size = 7), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 7), 
        legend.text = element_text(size = 5),
        legend.position = "bottom") 

xAxisBoxPlot_bm <- ggplot(data = BacSp_data.scores_mothers, aes(x = Timepoint, y = NMDS1)) +
  geom_boxplot(aes(fill=Timepoint), alpha=0.7, show.legend = FALSE) + 
  ylim(-1,1)+
  ggtitle("Shift in mother bacteriome composition over time") +
  scale_fill_manual(values=Peru1) +
  labs(tag = 'e') +
  coord_flip() +
  theme_void() +
  theme(plot.title = element_text(hjust=0.0, size=7), plot.tag = element_text(face="bold", vjust=-7, size=7))

Fig1E <- xAxisBoxPlot_bm / pbm + 
  plot_layout(heights = c(2,8)) + guides(color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
                                         alpha="none")

Fig1BC <- (Fig1B | Fig1C) + plot_layout(ncol=2, guides = "collect") & 
  theme(legend.position = "bottom", 
        legend.box.margin = margin(-15,0,0,0))

Fig1DE <- (Fig1D | Fig1E) + plot_layout(ncol=2, guides = "collect") & 
  theme(legend.position = "bottom", 
        legend.box.margin = margin(-15,0,0,0))

Fig1BCDE <- plot_spacer() / (Fig1BC) / (Fig1DE) + plot_layout(ncol=1, heights = c(4,3,3)) & theme(legend.position = "bottom", 
                                                                                                         #axis.title = element_text(size=7),
                                                                                                         #axis.text = element_text(size=5), 
                                                                                                         plot.title = element_text(size=7),
                                                                                                         legend.title = element_text(size=7),
                                                                                                         legend.text = element_text(size=5), 
                                                                                                         plot.tag = element_text(size=7))

pdf('Figures/Figure1/Fig1BCDE_1412.pdf', width=18/2.54, height=21/2.54)
Fig1BCDE
dev.off()


##### Figure 2A #######
Fig2A <- ggplot(VLP_metadata, aes(Timepoint, viral_alpha_diversity, fill=Type)) + 
  labs (y="Shannon Diversity Index", x="Timepoint", tag='a',colour = "Type") + 
  ylim(0,10) + 
  geom_violin(aes(color=Type),alpha=0.1, show.legend = F) + 
  geom_sina(aes(colour=Type), size=0.3, show.legend = F) +
  ggnewscale::new_scale_color() + 
  geom_boxplot(aes(color=Type),width = 0.3, outlier.shape = NA, alpha=0.3, show.legend = F, lwd=0.3) + 
  scale_color_manual(values=c("#B31312", "#090580")) + 
  ggtitle("Alpha diversity of vOTUs over time") +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  geom_text(data=significance_annotation_F2A, x=3, y=10, aes(label=label), size=1.75) +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7), 
        plot.tag = element_text(face="bold", size=7, vjust=-4))

##### Figure 2B #######
Fig2B <- ggplot(MGS_metadata, aes(Timepoint, bacterial_alpha_diversity, fill=Type)) + 
  labs (y="Shannon Diversity Index", x="Timepoint", tag="b") + 
  ylim(0,5) + 
  geom_violin(aes(color=Type),alpha=0.1, show.legend = F) + 
  geom_sina(aes(colour=Type), size=0.3, show.legend = F) +
  ggnewscale::new_scale_color() + 
  geom_boxplot(aes(color=Type),width = 0.3, outlier.shape = NA, alpha=0.3, show.legend = F, lwd=0.3) + 
  scale_color_manual(values=c("#B31312", "#090580")) + 
  ggtitle("Alpha diversity of bacterial species over time") +
  theme_bw()+
  theme(legend.position="none") + 
  facet_grid(~Type, scales="free", space = "free") +
  geom_text(data=significance_annotation_F2B, x=3.5, y=4.25, aes(label=label), size=1.75) +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7),
        plot.tag = element_text(face="bold", size=7, vjust=-4))

##### Figure 2C #######
Fig2C <- ggplot(data=virstability_M1_stat[virstability_M1_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_sina(data=virstability_M1_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=virstability_M1_stat[virstability_M1_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=0.5)  + 
  labs(x="Timepoint", y="N detected vOTUs", fill='vOTUs', color="vOTUs", tag="c") + 
  ggtitle("Infant virome") +
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
  scale_colour_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from M1', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained from M1', 'All')) +
  guides(#color=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5),
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

##### Figure 2D #######
Fig2D <- ggplot(virstability_M1_abundance_stat[virstability_M1_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position="identity") +
  geom_sina(data=virstability_M1_abundance_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=virstability_M1_abundance_stat[virstability_M1_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap, 
                    colour=Condition), 
                width=0.3, alpha=0.9, size=0.5)  +
  labs(x="Timepoint", y="Relative abundance",fill='vOTUs', color='vOTUs', tag = "d") + 
  ggtitle("Infant virome") +
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
  scale_colour_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from M1', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained from M1', 'All')) +
  guides(#color="none", 
    fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
    alpha="none")

##### Figure 2E #######
Fig2E <- ggplot(virstability_P7_stat[virstability_P7_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_sina(data=virstability_P7_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=virstability_P7_stat[virstability_P7_stat$Condition!='Not_retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap,
                    colour=Condition), 
                width=0.3, alpha=0.9, size=0.5)  +
  labs(x="Timepoint", y="N detected vOTUs",fill='vOTUs', color='vOTUs', tag="e") + 
  ggtitle("Mother virome") +
  theme_bw() + 
  theme(axis.title=element_text(size=7),
        axis.text=element_text(size=5),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        plot.title = element_text(hjust=0.5, size=7),
        legend.key.size=unit(0.7, "line")) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P7', 'All')) + 
  scale_color_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from P7', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3)) +
  guides(#color="none", 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

##### Figure 2F #######
Fig2F <- ggplot(virstability_P7_abundance_stat[virstability_P7_abundance_stat$Condition!='Not_retained',], aes(Timepoint, mean_bootstrap, fill=Condition, alpha=Condition)) + 
  geom_bar(color='grey', stat = "identity", position = "identity") +
  geom_sina(data=virstability_P7_abundance_dots, aes(Timepoint, value, color=Condition, group=Timepoint), size=0.3) +
  geom_errorbar(data=virstability_P7_abundance_stat[virstability_P7_abundance_stat$Condition=='Retained',], 
                aes(x=Timepoint, 
                    ymin=q1_bootstrap, 
                    ymax=q2_bootstrap,
                    colour=Condition), 
                width=0.3, alpha=0.9, size=0.5)  +
  labs(x="Timepoint",y="Relative abundance",fill='vOTUs', color='vOTUs', tag = "f") + 
  ggtitle("Mother virome") +
  theme_bw() + 
  theme(axis.title=element_text(size=7),
        axis.text=element_text(size=5),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5), 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        plot.title = element_text(hjust=0.5, size=7),
        legend.key.size=unit(0.7, "line")) + 
  scale_fill_manual(values = c(Demuth[5:6], "#FFFFFF"), labels=c('Retained from P7', 'All')) + 
  scale_color_manual(values = c(Demuth[c(3,8)]), labels=c('Retained from P7', 'All')) + 
  scale_alpha_manual(values=c(0.9,0.3), labels=c('Retained', 'Total space')) +
  guides(#color="none", 
         fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         alpha="none")

##### Figure 2G #######
Fig2G <- ggplot(plot_vir_frac_infants_per_sample_N, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Infant samples",y="% of sample composition", fill='Fraction', tag="g") +
  ggtitle("Infant virome") +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPV', 'TDV')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        plot.title = element_text(hjust=0.5, size=7),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5))

##### Figure 2H #######
Fig2H <- ggplot(p_vir_frac_infants_per_sample_ab, aes(Sample, value, fill=variable)) + 
  geom_bar(position = "fill", stat = "identity", width=1) + 
  labs(x="Infant samples",y="% relative abundance",fill='Fraction', tag="h") +
  ggtitle("Infant virome") +
  facet_grid(~Infant, space = "free_x", scales = "free_x") +
  theme_bw() + 
  scale_fill_manual(values=Veronese[5:4], labels=c('PPV', 'TDV')) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.y=element_text(size=5), 
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        axis.title =element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.position = "bottom", 
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        panel.spacing = unit(0.1, "lines"),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        plot.title = element_text(hjust=0.5, size=7),
        legend.key.size=unit(0.7, "line")) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE, title.position = 'top', title.hjust = 0.5))

##### Figure 2I #######

Fig2I <- ggplot(HostGen_average, aes(variable, value, fill=taxa, alpha=Sequencing)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Type, scales = 'free_x', nrow = 1) + 
  labs(x="Timepoint",y="Rescaled relative abundance",fill='Bacterial host genera', tag="i", alpha="Virome\nsource") +
  ggtitle("vOTUs aggregated at host genus level") +
  theme_bw() + 
  scale_fill_manual(values=Gen_palette[Gen_palette$Genus %in% unique(HostGen_average$taxa), ]$Color) +
  scale_alpha_manual(values = c(1,0.6), labels=c('VLP metavirome', 'MGS metavirome')) +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7),
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=4, byrow=F, title.position = 'top', title.hjust = 0.5,label.theme = element_text(face = "italic", size=5)), 
         alpha=guide_legend(nrow=2, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

##### Figure 2J #######
Fig2J <- ggplot(BacGen_average, aes(variable, value, fill=taxa)) + 
  geom_bar(position = "fill", stat = "identity") + 
  facet_wrap(~Type, scales = 'free_x', nrow = 1) + 
  xlab(label = "Timepoint") + 
  ylab(label = "Rescaled relative abundance") + 
  labs(fill='Bacterial host genera', tag="j") +
  ggtitle("Bacterial genera") +
  theme_bw() + 
  scale_fill_manual(values=Gen_palette[Gen_palette$Genus %in% unique(BacGen_average$taxa), ]$Color) +
  scale_alpha_manual(values = c(0.8, 1), labels=c('VLP metavirome', 'MGS metavirome')) +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7),
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill=guide_legend(nrow=4, byrow=F, title.position = 'top', title.hjust = 0.5,  label.theme = element_text(face = "italic", size=5)), 
         alpha=guide_legend(nrow=2, byrow=TRUE, title.position = 'top', title.hjust = 0.5, label.theme = element_text(size=5)))

Fig2AB <- (Fig2A | Fig2B ) + plot_layout(ncol=2,guides="keep") & theme(legend.position = "left")

Fig2CDEF <- (Fig2C | Fig2D | Fig2E | Fig2F) + plot_layout(ncol=4,guides="collect") & theme(legend.position = "bottom")

Fig2GH <- (Fig2G | Fig2H) +  plot_layout(ncol=2,guides="collect") & theme(legend.position = "bottom")

Fig2IJ <- (Fig2I | Fig2J ) + plot_layout(ncol=2,guides="collect") & theme(legend.position = "bottom")

Fig2 <- Fig2AB / Fig2CDEF / Fig2GH / Fig2IJ  + plot_layout(nrow=4, guides = "collect") & theme(#plot.title = element_text(size=7),
                                                                                                #                    legend.title = element_text(size=7),
                                                                                                 #                   legend.text = element_text(size=5),
                                                                                                                    legend.margin=margin(0,0,0,0),
                                                                                                                    legend.box.margin=margin(-8,-10,-6,-10))


pdf('Figures/Figure2/Fig2_1812.pdf', width=18/2.54, height=21/2.54)
Fig2
dev.off()

##### Figure 3 #######

##### Figure 3A #######
Fig3A <- ggplot(VLP_metadata, aes(Timepoint, temperate_RA, fill=Type)) + 
  ylim(0, 103) +
  labs (y="Relative abundance of\nactive temperate phages", x="Timepoint", tag="a") + 
  geom_boxplot(aes(color=Type),width=0.5,outlier.shape = NA, alpha=0.4, show.legend = F, lwd=0.3) + 
  geom_sina(aes(color=Type), size=0.6, show.legend = F) +
  facet_grid(~Type, scales = 'free_x') +
  geom_text(data=significance_annotation_F3A, x=3, y=102, aes(label=label), size=1.75) +
  theme_bw()+
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        strip.text.x = element_text(size = 6),
        strip.background = element_rect(color="black", fill=NA),
        legend.text = element_text(size=7),
        legend.title = element_text(size=5),
        plot.title = element_text(hjust=0.5, size=7), 
        plot.tag = element_text(face="bold", size=7))

##### Figure 3B #######
Fig3B <- ggplot(VLP_metadata[!is.na(VLP_metadata$infant_ever_never_breastfed),], aes(Timepoint, viral_alpha_diversity, fill=infant_ever_never_breastfed)) + 
  labs (y="vOTUs Shannon Diversity Index", x="Timepoint", tag="b") + 
  geom_boxplot(aes(color=infant_ever_never_breastfed), outlier.shape = NA, alpha=0.5, lwd=0.3) + 
  geom_sina(aes(color=infant_ever_never_breastfed), size=0.6) +
  annotate(geom="text", x=3, y=6, label="p-value=0.02, beta=0.9", size=1.75) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.tag = element_text(face="bold", size=7),
        legend.key.size=unit(0.7, "line")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Breastfed', 'Exclusively\nformula fed'),
                    values=c("#EA5455", "#002B5B")) +
  scale_color_manual(name = "Feeding mode", 
                     labels=c('Breastfed', 'Exclusively\nformula fed'),
                     values=c("#EA5455", "#002B5B")) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         color=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

##### Figure 3C #######
Fig3C <- ggplot(VLP_metadata[!is.na(VLP_metadata$infant_ever_never_breastfed),], aes(Timepoint, temperate_richness, fill=infant_ever_never_breastfed)) + 
  labs (y=expression(Log["10"]~"N active temperate phages"), x="Timepoint", tag="c") + 
  geom_boxplot(aes(color=infant_ever_never_breastfed), outlier.shape = NA, alpha=0.5, lwd=0.3) + 
  geom_sina(aes(color=infant_ever_never_breastfed), size=0.6) +
  annotate(geom="text", x=3, y=250, label="p-value=0.01, beta=15.0", size=1.75) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.tag = element_text(face="bold", size=7),
        legend.key.size=unit(0.7, "line")) +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('Breastfed', 'Exclusively\nformula fed'),
                    values=c("#EA5455", "#002B5B")) +
  scale_color_manual(name = "Feeding mode", 
                     labels=c('Breastfed', 'Exclusively\nformula fed'),
                     values=c("#EA5455", "#002B5B")) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         color=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

##### Figure 3D #######
lab_italics <- paste0("italic('", rownames(diff_ab_temperate), "')")
selectLab_italics = paste0(
  "italic('",
  c('Bacteroides fragilis','Phocaeicola vulgatus','Bacteroides caccae', 'Bacteroides salyersiae'),
  "')")
Fig3D <- EnhancedVolcano(diff_ab_temperate,
                         lab = lab_italics, parseLabels = TRUE,
                         x = 'log2FC',
                         y = 'FDR', 
                         ylim = c(0, -log10(1e-4)),
                         xlim = c(-6,6),
                         title = '',
                         pCutoff = 0.05,
                         #FCcutoff = 1.75,
                         drawConnectors = TRUE,
                         labSize = 1.75,
                         widthConnectors = 0.75,
                         boxedLabels = TRUE, 
                         colAlpha = 0.9,
                         colGradient = c('#B31312', '#2B2A4C'),
                         caption = "",
                         captionLabSize = 0, 
                         titleLabSize = 1, 
                         subtitleLabSize = 0,
                         subtitle = "",
                         selectLab = selectLab_italics, #c('Bacteroides fragilis', 'Phocaeicola vulgatus','Bacteroides caccae', 'Bacteroides salyersiae'),
                         axisLabSize = 5) + 
  ggplot2::labs(y=expression(-Log["10"]~FDR), x=expression(Log["2"]~"fold change of N temperate phages per bacterial host"), tag="d", color="FDR") +
  ggplot2::annotate(geom="text", x=-3.7, y=4, label="Enriched in breastfed", size=1.75) +
  ggplot2::annotate(geom="text", x=3.7, y=4, label="Enriched in formula fed", size=1.75) +
  ggplot2::ggtitle("Differential active temperate phage counts\nper bacterial host")+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.tag = element_text(face="bold", size=7, vjust = -4),
        plot.title = element_text(size=7, hjust=0.5, face='plain'),
        legend.position = "bottom",
        legend.key.size=unit(0.7, "line"))

##### Figure 3E #######

Fig3E <- ggplot(sharedness_timepoints_melt, aes(Timepoint, value*100, fill=variable)) +
  geom_boxplot(aes(color=variable), outlier.shape = NA, alpha=0.5, lwd=0.3) +
  labs (y="% infant vOTUs", x="Infant Timepoint", tag='e') + 
  geom_sina(aes(color=variable), size=0.6,alpha=0.5) +
  annotate(geom="text", x=3, y=85, label="p-value=0.04, beta=3.3e-02", size=1.75) +
  ylim(0,88) +
  ggtitle("Shared with active maternal virome\n") +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.title = element_text(size=7, hjust=0.5),
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.key.size=unit(0.7, "line")) +
  scale_fill_manual(name = "Maternal Timepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#F6830F", "#0E918C")) +
  scale_color_manual(name = "Maternal Timepoints", 
                     labels=c('Pre birth', 'Post birth'),
                     values=c("#F6830F", "#0E918C")) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         color=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

##### Figure 3F #######
Fig3F <- ggplot(sharedness_timepoints_iVM_mVM_melt, aes(Timepoint, value*100, fill=variable)) +
  geom_boxplot(aes(color=variable), outlier.shape = NA, alpha=0.5, lwd=0.3) +
  labs (y="% infant prophage-inclusive vOTUs", x="Infant Timepoint", tag="f") + 
  geom_sina(aes(color=variable), size=0.6,alpha=0.5) +
  annotate(geom="text", x=3, y=85, label="p-value=0.001, mean increase 4.9%", size=1.75) +
  ylim(0,88) +
  ggtitle("Shared with prophage-inclusive\nmaternal virome") +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),plot.title = element_text(size=7, hjust=0.5),
        plot.tag = element_text(face="bold", size=7, vjust=-4),
        legend.key.size=unit(0.7, "line")) +
  scale_fill_manual(name = "Maternal Timepoints", 
                    labels=c('Pre birth', 'Post birth'),
                    values=c("#F6830F", "#0E918C")) +
  scale_color_manual(name = "Maternal Timepoints", 
                     labels=c('Pre birth', 'Post birth'),
                     values=c("#F6830F", "#0E918C"))+ 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5), 
         color=guide_legend(nrow=1, byrow=TRUE, title.position = 'top', title.hjust = 0.5))

Fig3BC <- (Fig3B | Fig3C ) + plot_layout(ncol=2,guides="collect") & theme(legend.position = "bottom")

Fig3ABC <- (Fig3A | Fig3BC) + plot_layout(widths=c(4.5,5.5),guides="keep") & theme(legend.position = "bottom")

Fig3EF <- (Fig3E | Fig3F ) + plot_layout(ncol=2,guides="collect") & theme(legend.position = "bottom")

Fig3DEF <- (Fig3D | Fig3EF) + plot_layout(widths=c(4, 6), guides="keep") & theme(legend.position = "bottom")

Fig3 <- (Fig3ABC / Fig3DEF) + plot_layout(guides="keep") & theme(plot.title = element_text(size=7),
                                                                 legend.position = "bottom",
                                                                   legend.title = element_text(size=7),
                                                                   legend.text = element_text(size=5),
                                                                   legend.margin=margin(-5,0,0,0),
                                                                   legend.box.margin=margin(-8,-10,-6,-10))

pdf('Figures/Figure3/Fig3_1812.pdf', width=18/2.54, height=14/2.54)
Fig3 
dev.off()

##### Figure 4A #######
pvd_0 <- ggplot(plot_distances_select, aes(vector4analysis, easy_name, fill=factor4analysis)) + 
  labs (y="Viruses", x=expression(Log["10"]~"Kimura distance"), fill="Kinship", color="Kinship", tag='a') + 
  geom_sina(aes(color=factor4analysis), size=0.2, alpha=0.8) +
  geom_boxplot(outlier.shape = NA, alpha=0.3, lwd=0.3) +
  theme_bw()+
  scale_x_log10() +
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        #strip.text.x = element_text(size = 12),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        #legend.position = "none",
        plot.tag=element_text(face='bold', size=7),
        legend.key.size=unit(0.7, "line"),
        plot.margin = unit(c(0,-1,0,0), "cm")) +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  scale_color_manual(labels = c("Related", "Unrelated"), 
                     values=c("#17B971", "#7d8bb1")) +
  geom_text(data=plot_distances_select[plot_distances_select$factor4analysis=='Related',], aes(label = significance_level, x = Inf, vjust=2, y = easy_name), size = 2, angle=270) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5), 
         color=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5))

pvd_1 <- ggplot(N_pairwise_distance, aes(value, ContigID_easy, fill=variable) ) + 
  geom_bar(stat='identity', position='dodge', color="black", alpha=0.5, lwd=0.3) +
  scale_x_log10(breaks=c(10^0, 10^1, 10^3)) + 
  labs (y="", x=expression(Log["10"]~"N distances"),fill="Kinship") +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        #strip.text.x = element_text(size = 12),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,-0.3), "cm")) +
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  guides(fill = "none") 

Fig4A <- pvd_0  + pvd_1 +
  plot_layout(ncol = 2, nrow = 1, guides="collect", widths = c(8.5,1.5)) + 
  plot_annotation(title = "") & theme(legend.position = "bottom",
                                      plot.title = element_text(size=0),
                                      legend.title = element_text(size=7),
                                      legend.text = element_text(size=5),
                                      legend.margin=margin(0,0,0,0),
                                      legend.box.margin=margin(1,-5,-5,-5)) & guides(color="none")

##### Figure 4B #######
Fig4B <- ggplot(for_plot, aes(value*100, ContigID_easy, fill=variable)) + 
  geom_col(position = 'dodge', alpha=1) +
  labs (y="Viruses", x="% mother-infant sample pairs sharing virus", tag="b") + 
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7), 
        legend.position = 'bottom',
        plot.tag = element_text(face='bold', size=7),
        legend.key.size=unit(0.7, "line")) +
  labs(fill="") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  geom_text(data=for_plot[for_plot$variable=='Perc_related_pairs_transmitted',], aes(label = significance_level_transmission, x = 105, y = ContigID_easy), size = 1.75, angle=270) + 
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5))

##### Figure 4C #######
plot_distances_select_bac$easy_species_names <- paste0('*', gsub('_.*', '', plot_distances_select_bac$species_names), '*', '<br>', '(', plot_distances_select_bac$bacterium_name, ')')

pbd_0 <- ggplot(plot_distances_select_bac, aes(vector4analysis, ord, fill=factor4analysis)) + 
  labs (y="Bacterial strains", x=expression(Log["10"]~"Kimura distance"), tag="c") + 
  geom_sina(aes(color=factor4analysis), size=0.2, alpha=0.8) +
  geom_boxplot(outlier.shape = NA,alpha=0.3, lwd=0.3) +
  theme_bw()+
  scale_y_discrete(labels = setNames(as.character(plot_distances_select_bac$easy_species_names), plot_distances_select_bac$ord)) +
  scale_x_log10(breaks=c(1e-01, 1e+00, 1e+01, 1e+02), labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text=element_text(size=5),
        axis.text.y = ggtext::element_markdown(),
        axis.title=element_text(size=7),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        plot.tag = element_text(face="bold", size=7), 
        legend.key.size=unit(0.7, "line"),
        plot.margin = unit(c(0,-1,0,0), "cm")) +
  labs(fill="Kinship", color="") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) + 
  scale_color_manual(labels = c("Related", "Unrelated"), 
                     values=c("#17B971", "#7d8bb1")) +
  geom_text(data=plot_distances_select_bac[plot_distances_select_bac$factor4analysis=='Related',], aes(label = significance_level, x = Inf, vjust=2, y = ord), size = 1.75, angle=270) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5), 
         color="none")

pbd_1 <-  ggplot(N_pairwise_distance_bac, aes(value, ord, fill=variable) ) + 
  geom_bar(stat='identity', position='dodge', color="black", alpha=0.5, lwd=0.3) +
  scale_y_discrete(labels = setNames(as.character(N_pairwise_distance_bac$species_name), N_pairwise_distance_bac$ord)) +
  scale_x_log10(breaks=c(10^0, 10^1, 10^3)) + 
  labs (y="", x=expression(Log["10"]~"N distances")) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        #strip.text.x = element_text(size = 12),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,-0.3), "cm")) +
  labs(fill="")+
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  guides(fill = "none")

Fig4C <- pbd_0 + pbd_1 +
  plot_layout(ncol = 2, nrow = 1, guides="collect", widths = c(8, 2)) + 
  plot_annotation(title = "") & theme(legend.position = "bottom",
                                      plot.title = element_text(size=0),
                                      legend.title = element_text(size=7),
                                      legend.text = element_text(size=5),
                                      legend.margin=margin(0,0,0,0),
                                      legend.box.margin=margin(1,-5,-5,-5)) & guides(color="none")

##### Figure 4D #######
for_plot_bac$easy_species_names <- paste0('*', gsub('_.*', '', for_plot_bac$species_names), '*', '<br>', '(', for_plot_bac$Host_SGB, ')')

Fig4D <- ggplot(for_plot_bac, aes(value*100, ord, fill=variable)) + 
  geom_col(position = 'dodge', alpha=0.8) +
  labs (y="Bacterial strains", x="% mother-infant sample pairs sharing strain", tag="d") + 
  scale_y_discrete(labels = setNames(as.character(for_plot_bac$easy_species_names), for_plot_bac$ord)) +
  theme_bw()+
  theme(axis.text=element_text(size=5), 
        axis.title=element_text(size=7),
        axis.title.x = element_text(size=7, hjust=0.7),
        axis.text.y = ggtext::element_markdown(),
        legend.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.position = "bottom",
        plot.tag = element_text(face="bold", size=7),
        legend.key.size=unit(0.7, "line")) +
  labs(fill="") + 
  geom_stripes(odd = "#33333333", even = "#00000000") +
  scale_fill_manual(labels = c("Related", "Unrelated"), 
                    values=c("#17B971", "#7d8bb1")) +
  geom_text(aes(label = significance_level_transmission, x = 105, y = ord), size = 1.75, angle=270) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE, title.position = 'left', title.hjust = 0.5))

Fig4AB <- (pvd_0 | pvd_1 | Fig4B ) + plot_layout(ncol=3, widths = c(4,2,4),guides="collect") & theme(legend.position = "bottom", 
                                                                                                  plot.margin = margin(0,2,0,0),
                                                                                                  legend.margin=margin(0,0,0,0),
                                                                                                  legend.box.margin=margin(-10,0,0,0)) & guides(color="none")
Fig4CD <- (pbd_0 | pbd_1 | Fig4D) + plot_layout(ncol=3, widths = c(4,2.2,3.8),guides="collect") & theme(legend.position = "bottom", 
                                                                                                     plot.margin = margin(-9,2,-2.5,0),
                                                                                                     legend.margin=margin(0,0,0,0),
                                                                                                     legend.box.margin=margin(-10,0,0,0)) & guides(color="none")

Fig4ABCD <- wrap_elements(Fig4AB) / wrap_elements(Fig4CD)+ plot_layout(guides="keep", widths = c(4,4), heights = c(4.95, 5.05)) & theme(legend.position = "bottom",
                                                                                                               axis.text = element_text(size=7),                                       
                                                                                                               plot.title = element_text(size=0),
                                                                                                               legend.title = element_text(size=7),
                                                                                                               legend.text = element_text(size=5))

pdf('Figures/Figure4/Fig4_1512.pdf', width=18/2.54, height=21/2.54)
Fig4ABCD
dev.off()

##### Figure 5 #######

##### Figure 5A #####
pdf('./Figures/Figure5/Fig5A.pdf', width=21/2.54, height=19/2.54)
corrplot(bgcolors, na.label = "square", addgrid.col = NA, na.label.col = "#FC3C3C", tl.col = "white")
corrplot(co_transmission_cor,
         bg=NA,
         p.mat=co_transmission_FDR, 
         cl.lim=c(0,1), is.corr = F,
         insig = "blank",
         na.label = "X", 
         na.label.col = "lightgrey", 
         tl.col='black', 
         add = T,
        #outline=F,
         addgrid.col="darkgrey",
         number.cex=1, #X signs
         cl.cex=1, # X-axis text
         tl.cex = 1, # Y-axis text
         col=colorRampPalette(c("#FFFFFF","#828AB1","#041562"))(100), cl.pos="b")
dev.off()

##### Figure 5B #####
Fig5B <- ggplot(cotransmission_stat, aes(x=Var2, y=prop, width=paired.count, fill=Var1)) + 
  labs(x="Virus - predicted host pair", y="", fill="Transmission\npattern", tag="b") + 
  geom_bar(stat = "identity", position = "fill", colour = "black") + 
  geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5), size=2) +
  facet_grid(~Var2, scales = "free_x", space = "free_x")  +
  scale_fill_manual(labels=c("Co-transmitted", "Not co-transmitted"),values=c(Kandinsky[1],Kandinsky[2])) + 
  scale_y_continuous(breaks = c(0.28, 0.78), labels = c('Co-transmitted', 'Not co-transmitted')) +
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text=element_text(size=5), 
        axis.text.y = element_text(size=5, angle=90, hjust=0.5),
        axis.title=element_text(size=7),
        panel.spacing.x = unit(0, "npc"),
        legend.text = element_text(size=5, hjust=0.5),
        legend.title = element_text(size=7),
        plot.tag = element_text(face="bold", size=7),
        legend.position = "right",
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0)) + 
  guides(fill = guide_legend(nrow=2, title.position = "top"))

pdf("Figures/Figure5/Fig5B.pdf", width=8/2.54, height=7.5/2.54)
Fig5B
dev.off()
##### Figure 5C #####
ptv1 <- phylo_L_85266_LS0 + 
  geom_tippoint(aes(shape = Type, fill = FAMILY, color=Sequencing), size = 1.75, stroke=0.5) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L_85266_LS0_label_colors_family$color, L_85266_LS0_label_colors_family$L_85266_LS0_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.3, size=1.75) + 
  scale_color_manual(values=setNames(L_85266_LS0_label_colors$color, L_85266_LS0_label_colors$L_85266_LS0_SAMs)) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(aes(shape = Type, fill = FAMILY, colour = Sequencing), size = 1.75, stroke=0.5) +
  scale_colour_manual(values = c("black", "grey"), labels=c("VLP metavirome", "MGS metavirome\nTotal metagenome")) + 
  geom_treescale() +
  guides(#shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21,  stroke = NA))) +
  ggtitle("L85266_LS0") +
  xlim(0, 8) +
  labs(tag="c", colour="Strain source", shape="Family member", fill="Family") +
  theme(plot.tag = element_text(face="bold", size=7),
        plot.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0),
        plot.margin = margin(0,0,0,0))

ptb1 <- phylo_Host_L_85266_LS0 +
  geom_tippoint(aes(shape = Type, fill = FAMILY, color = Sequencing), size = 1.75, stroke=0.5) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L_85266_LS0_label_colors_family$color, L_85266_LS0_label_colors_family$L_85266_LS0_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.3, size=1.75) + 
  scale_color_manual(values=setNames(Host_L_85266_LS0_label_colors_bac$color, Host_L_85266_LS0_label_colors_bac$Host_L_85266_LS0)) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(aes(shape = Type, fill = FAMILY, colour = Sequencing), size = 1.75, stroke=0.5) +
  scale_colour_manual(values = c("grey","black"), labels=c("MGS metavirome\nTotal metagenome","VLP metavirome")) + 
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, stroke = NA))) +
  ggtitle("*Bacteroides uniformis*") + 
  xlim(0, 4) +
  labs(colour="Strain source", shape="Family member", fill="Family") +
  theme(plot.title = ggtext::element_markdown(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0),
        plot.margin = margin(0,0,0,0))


Fig5C <- (ptv1 | ptb1 ) + 
  plot_layout(ncol=2, guides="collect") & 
  theme(legend.position = "right") & 
  guides(shape = guide_legend(override.aes = list(color = NULL), nrow=3, title.position = 'top', title.hjust = 0.5),
         fill = guide_legend(override.aes = list(shape = 21, stroke=NA), nrow=7, title.position = 'top', title.hjust = 0.5),
         color = guide_legend(override.aes = list(shape=16), nrow=2, title.position = 'top', title.hjust = 0.5))


pdf('Figures/Figure5/Fig5C.pdf', width=18/2.54, height=13.5/2.54)
Fig5C 
dev.off()

##### Figure 6 #######

##### Figure 6A #######
Fig6A_v <- phylo_L34922_LS1 + 
  geom_tippoint(aes(shape = Type, fill = FAMILY, color=Sequencing), size = 1.75, stroke=0.5) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L34922_LS1_label_colors_family$color, L34922_LS1_label_colors_family$L34922_LS1_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.9, size=1.75) + 
  scale_color_manual(values=setNames(L34922_LS1_label_colors$color, L34922_LS1_label_colors$L34922_LS1_SAMs)) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(aes(shape = Type, fill = FAMILY, colour = Sequencing), size = 1.75, stroke=0.5) +
  scale_colour_manual(values = c("grey","black"), labels=c("MGS metavirome\nTotal metagenome","VLP metavirome")) + 
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, stroke = NA))) +
  ggtitle("L34922_LS1") +
  xlim(0, 20) +
  labs(tag="a", colour="Strain source", shape="Family member", fill="Family") +
  theme(plot.tag = element_text(face="bold", size=7),
        plot.title = element_text(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))

Fig6A_b <- phylo_Host_L34922_LS1 + 
  geom_tippoint(aes(shape = Type, fill = FAMILY, color = Sequencing), size = 1.75, stroke=0.5) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L34922_LS1_label_colors_family$color, L34922_LS1_label_colors_family$L34922_LS1_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.02, size=1.75) + 
  scale_color_manual(values=setNames(Host_L34922_LS1_label_colors_bac$color, Host_L34922_LS1_label_colors_bac$Host_L34922_LS1)) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(aes(shape = Type, fill = FAMILY, colour = Sequencing), size = 1.75, stroke=0.5) +
  scale_colour_manual(values = c("grey","black"), labels=c("MGS metavirome\nTotal metagenome","VLP metavirome")) + 
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, stroke = NA))) +
  ggtitle("*Bifidobacterium bifidum* <br> (SGB17256)") + 
  labs(colour="Strain source", shape="Family member", fill="Family") +
  xlim(0, 0.8) +
  theme(plot.title = ggtext::element_markdown(size=7),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        legend.key.size=unit(0.7, "line"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))

Fig6A <- (Fig6A_v | Fig6A_b) + 
  plot_layout(ncol=2, guides="collect") & 
  theme(legend.position = "right") & 
  guides(shape = guide_legend(override.aes = list(color = NULL), nrow=3, title.position = 'top', title.hjust = 0.5),
         fill = guide_legend(override.aes = list(shape = 21, stroke=NA), nrow=7, title.position = 'top', title.hjust = 0.5),
         color = guide_legend(override.aes = list(shape=16), nrow=2, title.position = 'top', title.hjust = 0.5))

##### Figure 6B #######
Fig6B <- ggplot(Host_L34922_LS1_copresence) +
  geom_polygon(aes(x = x, y = y, fill = MGS, group = interaction(Type, Timepoint)), color = "black", show.legend = F) +
  scale_fill_manual(values=c('white', '#D72323', '#7F7F7F'),  limits = c("0", "1", "NA")) +
  ggnewscale::new_scale_fill() +
  labs(x="Timepoint", y="", tag="b") +
  geom_polygon(aes(x = xdown, y = ydown, fill = VLP, group = interaction(Type, Timepoint)), color = "black") +
  scale_fill_manual(values=c('white', '#D72323', '#7F7F7F'),  limits = c("0", "1", "NA")) +
  ggtitle("Presence patterns<br>of L34922_LS1 and<br>*Bifidobacterium bifidum*") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                     labels = levels(Host_L34922_LS1_copresence$Timepoint)) +
  scale_y_continuous(breaks = c(1,3),
                     labels = c("Infant","Mother")) +
  coord_equal() + 
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=7),
        axis.text = element_text(size=5),
        legend.position = "right",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        plot.tag=element_text(face="bold", size=7),
        plot.title = ggtext::element_markdown(size=7),
        legend.key.size=unit(0.7, "line")) + 
  guides(fill= guide_legend("Presence", title.position = "top", title.hjust = 0.5))

##### Figure 6C #######
Fig6C <- plot_synteny(L34922_LS1_to_B.bifidum_patch_ali, 
                      q_chrom="L34922_LS1", 
                      t_chrom="B. bifidum patched\ngenome fragment", 
                      centre=TRUE,
                      xlab="Genome coordinate") + 
  labs(tag="c") +
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(size=7), 
        axis.text.x = element_text(size=5), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.tag = element_text(face="bold", size=7)) +
  annotate(geom="text", x = 750000, y=1.5, label="ANI > 99%\ncoverage=100%\ne-value=0", size=1.75) + 
  annotate(geom="text", x = 500000, y=0.70, label="B. bifidum patched genome fragment", size=1.75) + 
  annotate(geom="text", x=500000, y=2.3, label="L34922_LS1", size=1.75)

##### Figure 6D #######
Fig6D <-  gggenomes(genes=bifidophages_genes[bifidophages_genes$seq_id=='L34922_LS1',], 
                    seqs=bifidophages_seqs[bifidophages_seqs$seq_id=='L34922_LS1',]) + 
  geom_seq() +         # draw contig/chromosome lines
  geom_bin_label(size=1.75) +   # label each sequence 
  geom_gene(aes(fill=Color)) +         # draw genes as arrow
  geom_gene_label(aes(label=VOG_name_ed), size=1.75, vjust = 1) +
  labs(tag='d', fill='Genome coordinate, kb\nGene annotation source') +
  #ggtitle("Gene annotation") +
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


##### Figure 6E #######
Fig6E <- ggplot(L34922_LS1_MGS_reads, aes(w_center/1000, value, group=variable, color=color_factor, size=width)) +
  geom_line() + 
  xlim(0, 150) + 
  labs(x="*B. bifidum* patched genome fragment coordinate, kb", y="Log<sub>10</sub> total metagenome<br>mean read depth", color="Timepoint", tag='e') +
  scale_y_log10() + 
  scale_size_manual(values=c(0.5, 0.1), guide=FALSE) +
  scale_color_manual(values=L34922_LS1_read_align_colors$Color) + 
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
  geom_vline(data=L34922_LS1_prophage_region, aes(xintercept=coordinate/1000, color="Prophage region"), size=0.5, linetype="dashed") + 
  scale_color_manual(name="",
                     values=c(`Prophage region`="black"), 
                     guide = guide_legend(order = 3)) 

##### Figure 6F #######
Fig6F <- ggplot(L34922_LS1_VLP_reads, aes(w_center/1000, value, group=variable, color=color_factor, size=width)) +
  geom_line(show.legend = F) + 
  xlim(0, 150) + 
  labs(x="*B. bifidum* patched genome fragment coordinate, kb", y="Log<sub>10</sub> VLP metavirome<br> mean read depth", color="Timepoint", tag='f') +
  scale_y_log10() + 
  scale_size_manual(values=c(0.5, 0.1), guide=FALSE) +
  scale_color_manual(values=L34922_LS1_read_align_colors$Color) + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.tag = element_text(face="bold", size=7),
        axis.title.x = ggtext::element_markdown(size=7),
        axis.title.y = ggtext::element_markdown(size=7),
        axis.text = element_text(size=5),
        legend.title = element_text(size=7),
        legend.text = element_text(size=5)) + 
  guides(color=guide_legend(nrow=4, title.position="top", title.hjust = 0.5)) + 
  ggnewscale::new_scale_color() +
  geom_vline(data=L34922_LS1_prophage_region, aes(xintercept=coordinate/1000, color="Prophage region"), size=0.5, linetype="dashed") + 
  scale_color_manual(name="",
                     values=c(`Prophage region`="black"), 
                     guide = guide_legend(order = 3)) 


Fig6BC <- (Fig6B | Fig6C) + plot_layout(ncol = 2, guides="keep") & 
  theme(legend.position = "right") #& 
# guides(shape = guide_legend(override.aes = list(color = NULL), nrow=3, title.position = 'top', title.hjust = 0.5),
#        fill = guide_legend(override.aes = list(shape = 21, color = "black"), nrow=7, title.position = 'top', title.hjust = 0.5),
#        color = guide_legend(override.aes = list(shape=16), nrow=2, title.position = 'top', title.hjust = 0.5))

Fig6EF <- (Fig6E | Fig6F) + plot_layout(ncol = 2, guides="collect") & 
  theme(legend.position = "bottom")

pdf('Figures/Figure6/Fig6_1712_ed.pdf', width=18/2.54, height=21/2.54)
(Fig6A / Fig6BC / Fig6D / Fig6EF) + plot_layout(heights=c(6,2,4,5),guides="keep")
dev.off()

pdf('Figures/Figure6/Fig6C.pdf', width=26/2.54, height=10/2.54)
(Fig6C) 
dev.off()

pdf('Figures/Figure6/Fig6D.pdf', width=18/2.54, height=8/2.54)
(Fig6D) 
dev.off()
