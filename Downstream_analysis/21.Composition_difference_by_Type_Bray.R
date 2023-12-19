setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore the differences in bacteriome and virome 
# communitues between infant and maternal samples at the 
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

library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)
##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_with_phenos.txt', sep='\t', header=T)

MGS_metadata$Type <- factor(MGS_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)

microbiome <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_with_phenos.txt', sep='\t', header=T)

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header = T)
VLP_metadata$Type <- factor(VLP_metadata$Type, levels=c('Infant', 'Mother'), ordered=T)

RPKM_virome_aggr <- read.table('02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_raw.txt', sep='\t', header=T)
##############################
# ANALYSIS
##############################

#### NMDS for bacteriome

# preapring phenotypes
for_bacplot <- MGS_metadata
row.names(for_bacplot) <- MGS_metadata$Short_sample_ID_bact
for_bacplot <- for_bacplot[colnames(microbiome),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & BACTERIAL ALPHA DIVERSITY
for_bacplot <- for_bacplot[,c("Type", "Clean_reads", "DNA_CONC", "bacterial_alpha_diversity", "metaphlan_unknown_perc")]

ord <- metaMDS(t(microbiome), distance = "bray", k=2)
en = envfit(ord, for_bacplot, permutations = 999, na.rm = TRUE)
en$factors #r2 0.7416; p-value < 0.001
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Type <- c(gsub('Type', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs$Species <- c('N clean reads', 'DNA concentration', 'Alpha diversity (bacteria)', '% unknown in Metaphlan4')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_bacplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/Bacterial_spp_Bray_NMDS.pdf', width=10/2.54, height=10/2.54)
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2)+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', size = 2) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

ggMarginal(gg, type="densigram", groupFill=T)
dev.off()

### calculating the significance of separation:
ord2 <- metaMDS(t(microbiome), distance = "bray", k=1)

colnames(MGS_metadata)

infants_phenos_MGS <- mixed_models_taxa(MGS_metadata, 
                                        "Short_sample_ID_bact", 
                                        as.data.frame(scores(ord2, "sites")), 
                                        c( "Type"), "dont_consider_time")
#### FOR SUPPLEMENTARY TABLE ####
write.table(infants_phenos_MGS, '05.MANUSCRIPT/Supplementary_tables/MM_bacteriome_NMDS1_Type.txt', sep='\t', quote=F)

#### NMDS for virome

# preapring phenotypes
for_virplot <- VLP_metadata
row.names(for_virplot) <- VLP_metadata$Short_sample_ID
for_virplot <- for_virplot[colnames(RPKM_counts_VLP),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & VIRAL ALPHA DIVERSITY
for_virplot <- for_virplot[,c("Type", "Clean_reads", "DNA_CONC", "viral_alpha_diversity", "bacterial_contamination_perc_reads")]

#### yes, I use the same object names and don't make a function out of it
ord <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=2)
en = envfit(ord, for_virplot, permutations = 999, na.rm = TRUE)
en$factors #r2 0.33; p-value < 0.001
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Type <- c(gsub('Type', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs <- spp.scrs[spp.scrs$Species!='Clean_reads',] # because it was not significant in envfit
spp.scrs$Species <- c('DNA concentration', 'Alpha diversity (viruses)', '% bacterial contamination')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_virplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/Viral_vOTUs_Bray_NMDS.pdf', width=10/2.54, height=10/2.54)
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2)+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', alpha = 0.8, size = 2) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

ggMarginal(gg, type="densigram", groupFill=T)
dev.off()


### calculating the significance of separation:
ord2 <- metaMDS(t(RPKM_counts_VLP), distance = "bray", k=1)

colnames(VLP_metadata)

infants_phenos_VLP <- mixed_models_taxa(VLP_metadata, 
                                        "Short_sample_ID", 
                                        as.data.frame(scores(ord2, "sites")), 
                                        c( "Type"), "dont_consider_time")
#### FOR SUPPLEMENTARY TABLE ####
write.table(infants_phenos_VLP, '05.MANUSCRIPT/Supplementary_tables/MM_virome_NMDS1_Type.txt', sep='\t', quote=F)


#### NMDS for viruses aggregated based on bacterial hosts at the species-level

# preapring phenotypes
for_viragrplot <- VLP_metadata
row.names(for_viragrplot) <- VLP_metadata$Short_sample_ID
for_viragrplot <- for_viragrplot[colnames(RPKM_virome_aggr),]
#CHOOSING TECHNICAL PHENOTYPES, SAMPLE TYPE & VIRAL ALPHA DIVERSITY
for_viragrplot<- for_viragrplot[,c("Type", "Clean_reads", "DNA_CONC", "viral_alpha_diversity", "bacterial_contamination_perc_reads")]

#### yes, I use the same object names and don't make a function out of it
ord <- metaMDS(t(RPKM_virome_aggr), distance = "bray", k=2)
en = envfit(ord, for_viragrplot, permutations = 999, na.rm = TRUE)
en$factors #r2 0.37; p-value < 0.001
en$vectors

centroids <- as.data.frame(scores(en, "factors"))
centroids$Type <- c(gsub('Type', '', row.names(centroids)))

spp.scrs <- as.data.frame(scores(en, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
spp.scrs$Type <- 'Infant'
spp.scrs <- spp.scrs[spp.scrs$Species!='Clean_reads',] # because it was not significant in envfit
spp.scrs$Species <- c('DNA concentration', 'Alpha diversity (viruses)', '% bacterial contamination')

### for plotting:
data.scores = as.data.frame(scores(ord, "sites"))
data.scores <- merge(data.scores, for_viragrplot, by='row.names')
row.names(data.scores) <- data.scores$Row.names
data.scores$Row.names <- NULL

pdf('./04.PLOTS/Viruses_aggregated_Bray_NMDS.pdf', width=10/2.54, height=10/2.54)
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 2, alpha=0.8) + 
  geom_point(data=centroids, aes(fill=Type),shape=23, size=4, color='black', ) + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2)+
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_label_repel(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), color='black', alpha = 0.8, size = 2) +
  theme_bw()+
  theme(axis.text=element_text(size=10),
        axis.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30"),
        legend.position = "bottom") + 
  scale_color_manual(values=c("#F8766D", "#00BFC4"))

ggMarginal(gg, type="densigram", groupFill=T)
dev.off()


### calculating the significance of separation:
RPKM_virome_aggr <- RPKM_virome_aggr[,VLP_metadata$Short_sample_ID]
if (identical(colnames(RPKM_virome_aggr), VLP_metadata$Short_sample_ID)) {
  viruses_aggr_permanova <- adonis2(t(RPKM_virome_aggr) ~ Type +viral_alpha_diversity + Clean_reads + DNA_CONC + bacterial_contamination_perc_reads, 
                               data=VLP_metadata,  
                               permutations = 999,
                               method="bray",
                               by = "margin")
}


###### OUTPUT #####
write.table(bacteria_permanova, '03a.RESULTS/Bacterial_spp_bray_PERMANOVA.txt', sep='\t', quote=F)
write.table(viruses_permanova, '03a.RESULTS/Viral_vOTUs_bray_PERMANOVA.txt', sep='\t', quote=F)
write.table(viruses_aggr_permanova, '03a.RESULTS/Viruses_aggregated_bray_PERMANOVA.txt', sep='\t', quote=F)
