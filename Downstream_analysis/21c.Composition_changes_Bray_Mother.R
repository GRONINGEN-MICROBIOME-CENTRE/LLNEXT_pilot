setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the dynamics of bacteriome and virome 
# communitues within infant samples at the 
# species and vOTUs levels USING BRAY-CURTIS DISTANCE
#############################################################

##############################
# Functions
##############################
mixed_models_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list, consider_time) {
  df <- metadata
  row.names(df) <- df[,ID]
  df<-merge(df, CLR_transformed_data, by='row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Prevalent= c(colnames(CLR_transformed_data))
  #pheno_list= phenotypes
  
  Overall_result_phenos =tibble() 
  
  for (Bug in Prevalent){
    if (! Bug %in% colnames(df)){ next }
    #Prevalence = sum(as.numeric(as_vector(select(df, Bug)) > 0)) / dim(df)[1]
    # print (c(Bug, Prevalence))
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      if (consider_time=='time_as_covariate') {
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + (1|Individual_ID)"), collapse="" )) 
      } else { # else is mainly for associating entities with time alone
        Model0 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + (1|Individual_ID)"), collapse="" )) 
      }
      
      lmer(Model0, df_pheno) -> resultmodel0
      base_model=resultmodel0
      
      if (consider_time=='time_as_covariate') {
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + infant_mode_delivery + Age_months + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      } else { # else is mainly for associating entities with time alone
        Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + ",pheno2, "+ (1|Individual_ID)"), collapse="" ))
      }
      
      lmer(Model2, df_pheno, REML = F) -> resultmodel2
      M = "Mixed"
      as.data.frame(anova(resultmodel2, base_model))['resultmodel2','Pr(>Chisq)']->p_simp
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = p_simp, Model_choice = M, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p=as.data.frame(Overall_result_phenos)
  p <- p[! duplicated(paste0(p$Pheno, p$Bug)),]
  p$FDR<-p.adjust(p$P, method = "BH")
  
  return(p)
  
}
##############################
# Loading libraries
##############################
library(vegan)

library(ggplot2)
library(ggExtra)
library(ggrepel)
library(ggforce)
library(MetBrewer)

library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)
##############################
# Input data
##############################
VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T)
MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T)
common_columns <- c("Universal_fecal_ID", "SAMPLE_ID", "Old_ID", 
                    "NG_ID", "NEXT_ID", "Type", 
                    "Timepoint", "FAM_ID", "Short_sample_ID", 
                    "Short_sample_ID_bact", "Individual_ID",
                    "DNA_CONC", "Clean_reads", "mother_age_years","Timepoint_continuous")

metadata <- rbind(VLP_metadata[,common_columns], MGS_metadata[,common_columns])
metadata$source <- lapply(metadata$Short_sample_ID, function(x){
  ifelse(length(grep('V', x))!=0, "VLP", "MGS")
}   )
metadata$source <- as.character(metadata$source)

metadata <- metadata[metadata$Type=="Mother",]

metadata$Timepoint <- factor(metadata$Timepoint, levels=c('P3','P7','B',"M1", "M2", "M3"), ordered = T)

#row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID_bact

microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
microbiome <- microbiome[,metadata[!is.na(metadata$Short_sample_ID_bact),]$Short_sample_ID_bact]
microbiome <- microbiome[rowSums(microbiome) !=0,]

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_VLP <- RPKM_counts_VLP[,metadata[metadata$source=="VLP",]$Short_sample_ID]
RPKM_counts_VLP <- RPKM_counts_VLP[rowSums(RPKM_counts_VLP)!=0,]

RPKM_counts_MGS <- read.table('02.CLEAN_DATA/RPKM_counts_MGS.txt', sep='\t', header=T)
RPKM_counts_MGS <- RPKM_counts_MGS[,metadata[metadata$source=="MGS" & metadata$Timepoint=="P3",]$Short_sample_ID]

RPKM_counts_combo <- merge(RPKM_counts_VLP, RPKM_counts_MGS, by="row.names", all.x=T)
row.names(RPKM_counts_combo) <- RPKM_counts_combo$Row.names
RPKM_counts_combo$Row.names <- NULL
RPKM_counts_combo[is.na(RPKM_counts_combo)] <- 0
RPKM_counts_combo <- RPKM_counts_combo[rowSums(RPKM_counts_combo)>0,]

#row.names(VLP_metadata) <- VLP_metadata$Short_sample_ID

Peru1 <- met.brewer('Peru1')
##############################
# ANALYSIS
##############################
#### NMDS for bacteriome
# preapring phenotypes
for_bacplot <- metadata[metadata$source=="MGS",]
row.names(for_bacplot) <- for_bacplot$Short_sample_ID_bact
for_bacplot <- for_bacplot[colnames(microbiome),]

# calculating ordination with 2 dimensions
ord <- metaMDS(t(microbiome), distance = "bray", k=2)
# calculating significance of factors and vectors
en = envfit(ord, for_bacplot[,c("Timepoint_continuous", "Timepoint")], permutations = 999, na.rm = TRUE)
en$factors 
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Timepoint <- c(gsub('Timepoint', '', row.names(centroids)))
centroids$Timepoint <- factor(centroids$Timepoint, levels = c('P3','P7','B', "M1", "M2", "M3"), ordered = T)

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Mother'
spp.scrs$Species <- c('Time')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_bacplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL


pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS_Timepoint_color_Mother.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Timepoint),shape=23, size=4, color='black', ) + 
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  scale_fill_manual(values=Peru1) +
  scale_color_manual(values=Peru1) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") 
dev.off()

##### SAVE FOR PATCHING PANELS ####
write.table(data.scores, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS_scores_BacSp_mothers.txt", sep='\t', quote = F)
write.table(centroids, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS_centroids_BacSp_mothers.txt", sep='\t', quote = F)
write.table(spp.scrs, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS_vectors_BacSp_mothers.txt", sep='\t', quote = F)
##### SAVE FOR PATCHING PANELS ####

#### association to time

# allow only 1 dimension
ord2 <- metaMDS(t(microbiome), distance = "bray", k=1)

mothers_phenos_MGS <- mixed_models_taxa(for_bacplot, 
                                    "Short_sample_ID_bact", 
                                    as.data.frame(scores(ord2, "sites")), 
                                    c("Timepoint_continuous", "Timepoint"), "dont_consider_time")
#### FOR SUPPLEMENTARY TABLE ####
write.table(mothers_phenos_MGS, '05.MANUSCRIPT/Supplementary_tables/MM_mother_bacteriome_NMDS1_over_time.txt', sep='\t', quote=F)

data.scores2 = as.data.frame(scores(ord2, "sites"))
data.scores2 <- merge(data.scores2, for_bacplot, by='row.names')
row.names(data.scores2) <- data.scores2$Row.names
data.scores2$Row.names <- NULL

# posthoc:
TukeyHSD((aov(NMDS1 ~ Timepoint, data=data.scores2)))


pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS1_Timepoint_color_Mother.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores2, aes(x = Timepoint, y = NMDS1, fill=Timepoint)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  geom_sina(aes(color=Timepoint), size=0.6,alpha=0.7) +
  scale_fill_manual(values=Peru1) +
  scale_color_manual(values=Peru1) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") +
        labs(fill='Timepoint', color='Timepoint') +
  coord_flip()
dev.off()

##### SAVE FOR PATCHING PANELS ####
write.table(data.scores2, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1E_NMDS1_scores_BacSp_mothers.txt", sep='\t', quote = F)
##### SAVE FOR PATCHING PANELS ####

#### NMDS for virome
# preparing phenotypes
for_virplot <- metadata[metadata$Short_sample_ID %in% colnames(RPKM_counts_VLP),]
row.names(for_virplot) <- for_virplot$Short_sample_ID
for_virplot <- for_virplot[colnames(RPKM_counts_VLP),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & BACTERIAL ALPHA DIVERSITY

ord <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=2)
en = envfit(ord, for_virplot[,c("Timepoint_continuous", "Timepoint")], permutations = 999, na.rm = TRUE)
en$factors 
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Timepoint <- c(gsub('Timepoint', '', row.names(centroids)))
centroids$Timepoint <- factor(centroids$Timepoint, levels = c('P7','B', "M1", "M2", "M3"), ordered = T)

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Mother'
spp.scrs$Species <- c('Time')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_virplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/vOTUs_Bray_NMDS_Timepoint_color_Mother.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Timepoint)) + 
  stat_ellipse(geom = "polygon", alpha = 0.2, aes(group = Timepoint, fill=Timepoint), linetype = 2) +
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Timepoint),shape=23, size=4, color='black') + 
  scale_fill_manual(values=Peru1) +
  scale_color_manual(values=Peru1) +
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#444444") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2, alpha=0.7) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom") 
dev.off()

##### SAVE FOR PATCHING PANELS ####
write.table(data.scores, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS_scores_vOTUs_mothers.txt", sep='\t', quote = F)
write.table(centroids, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS_centroids_vOTUs_mothers.txt", sep='\t', quote = F)
write.table(spp.scrs, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS_vectors_vOTUs_mothers.txt", sep='\t', quote = F)
##### SAVE FOR PATCHING PANELS ####


#### infants phenotypes

# allow only 1 dimension
ord2 <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=1)

mother_phenos_VLP <- mixed_models_taxa(for_virplot, 
                                    "Short_sample_ID", 
                                    as.data.frame(scores(ord2, "sites")), 
                                    c("Timepoint", "Timepoint_continuous"), "dont_consider_time")
##### FOR SUPPLEMENTARY #####
write.table(mother_phenos_VLP, "05.MANUSCRIPT/Supplementary_tables/MM_NMDS1_maternal_virome_over_time.txt", sep='\t', row.names = F)


data.scores2 = as.data.frame(scores(ord2, "sites"))
data.scores2 <- merge(data.scores2, for_virplot, by='row.names')
row.names(data.scores2) <- data.scores2$Row.names
data.scores2$Row.names <- NULL

pdf('./04.PLOTS/vOTUs_Bray_NMDS1_Timepoint_color_Mother.pdf', width=10/2.54, height=10/2.54)
ggplot(data = data.scores2, aes(x = Timepoint, y = NMDS1, fill=Timepoint)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5) +
  geom_sina(aes(color=Timepoint), size=0.6,alpha=0.7) +
  scale_fill_manual(values=Peru1) +
  scale_color_manual(values=Peru1) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 8, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 8, colour = "grey30"),
        legend.position = "bottom")
dev.off()

##### SAVE FOR PATCHING PANELS ####
write.table(data.scores2, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig1F_NMDS1_scores_vOTUS_mothers.txt", sep='\t', quote = F)
##### SAVE FOR PATCHING PANELS ####
###### OUTPUT #####
