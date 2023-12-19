setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore changes in abundance of gut bacterial
# species, vOTUs, and aggregated vOTUs at the level of 
# bacterial species in the developing infant gut 
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

linear_model_taxa <- function(metadata, ID, CLR_transformed_data, pheno_list) {
  
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
    Bug2 = paste(c("`",Bug, "`"), collapse="")
    for ( pheno in pheno_list){
      pheno2 = paste(c("`",pheno, "`"), collapse="")
      df[is.na(df[colnames(df) == pheno]) == F, ID] -> To_keep
      df_pheno = filter(df, !!sym(ID) %in% To_keep )
      
      Model2 = as.formula(paste( c(Bug2,  " ~ DNA_CONC + Clean_reads + ",pheno2), collapse="" ))
      lm(Model2, df_pheno) -> resultmodel2
      
      as.data.frame(summary(resultmodel2)$coefficients)[grep(pheno, row.names(as.data.frame(summary(resultmodel2)$coefficients))),] -> Summ_simple
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(P = Summ_simple$`Pr(>|t|)`, Bug =Bug, Pheno=pheno, Model="simple") -> temp_output
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
library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)

library(reshape2)
library(ggplot2)
library(ggforce)

##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata <- MGS_metadata[MGS_metadata$Type=='Infant',]
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID_bact


species <- read.table('02.CLEAN_DATA/Microbiome_species_unfiltred.txt', sep='\t', header=T)
species <- species[,MGS_metadata$Short_sample_ID_bact]
# filtering for presence in more than 5% of microbiome infant samples:
species_filt <- species[(rowSums(species!=0) > 0.10*ncol(species)), ]
species_filt <- as.data.frame(t(species_filt))
# CLR-transformation
my_pseudocount_normal=min(species_filt[species_filt!=0])/2
species_filt_CLR<-decostand(species_filt, "clr", pseudocount=my_pseudocount_normal)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata <- VLP_metadata[VLP_metadata$Type=='Infant',]
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
row.names(VLP_metadata) <- VLP_metadata$Short_sample_ID

### vOTUs:
RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_VLP <- RPKM_counts_VLP[,VLP_metadata$Short_sample_ID]
# filtering for presence in more than 5% of virome infant samples:
RPKM_counts_VLP_filt <- RPKM_counts_VLP[( rowSums(RPKM_counts_VLP!=0) > 0.10*ncol(RPKM_counts_VLP) ), ]
RPKM_counts_VLP_filt <- as.data.frame(t(RPKM_counts_VLP_filt))
#CLR-transformation
my_pseudocount_normal=min(RPKM_counts_VLP_filt[RPKM_counts_VLP_filt!=0])/2
RPKM_counts_VLP_filt_CLR <- decostand(RPKM_counts_VLP_filt, "clr", pseudocount=my_pseudocount_normal)


### vOTUs aggregated at the species level
vOTUs_assigned_BacSp <- read.table('02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_species.txt', sep='\t', header = T)
vOTUs_assigned_BacSp <- vOTUs_assigned_BacSp[,VLP_metadata$Short_sample_ID]
# filtering for presence in more than 5% of virome infant samples:
vOTUs_assigned_BacSp_filt <- vOTUs_assigned_BacSp[(rowSums(vOTUs_assigned_BacSp!=0) > 0.10*ncol(vOTUs_assigned_BacSp)),]
vOTUs_assigned_BacSp_filt <- as.data.frame(t(vOTUs_assigned_BacSp_filt))
#CLR-transformation
my_pseudocount_normal=min(vOTUs_assigned_BacSp_filt[vOTUs_assigned_BacSp_filt!=0])/2
vOTUs_assigned_BacSp_filt_CLR <- decostand(vOTUs_assigned_BacSp_filt, "clr", pseudocount=my_pseudocount_normal)

##############################
# ANALYSIS
##############################
species_phenos <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", species_filt_CLR, c("infant_place_delivery", 
                                                                                              "infant_ffq_feeding_mode_complex",
                                                                                              "infant_ever_never_breastfed"), 'time_as_covariate')
#### FOR SUPPLEMENTARY TABLE ####
write.table(species_phenos, '05.MANUSCRIPT/Supplementary_tables/MM_BacSp_RA_phenos_over_time_infant.txt', sep='\t', quote=F)

crossectional_solid_food <- linear_model_taxa(MGS_metadata[MGS_metadata$Timepoint=='M6',], "Short_sample_ID_bact", species_filt_CLR, "infant_food_solid_intro_M6")


place_of_delivery_bacteria <- species_filt_CLR[,unique(species_phenos[species_phenos$FDR < 0.05 & species_phenos$Pheno %in% c('infant_place_delivery'),]$Bug), drop=F]
colnames(place_of_delivery_bacteria) <- gsub('.*s__', '', colnames(place_of_delivery_bacteria))
place_of_delivery_bacteria <- merge(place_of_delivery_bacteria, MGS_metadata[,c("Timepoint", 'infant_place_delivery', 'infant_mode_delivery', 'Age_months', 'DNA_CONC', 'Clean_reads', 'NEXT_ID')], by='row.names')

pdf('./04.PLOTS/Species_Akkermansia_muciniphila_infants_birthplace_no_CS.pdf', width=12/2.54, height=9/2.54)
ggplot(place_of_delivery_bacteria, aes(Timepoint, Akkermansia_muciniphila, fill=infant_place_delivery)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Place of delivery", 
                    labels=c('Home', 'Hospital'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()

ggplot(place_of_delivery_bacteria, aes(infant_place_delivery, Akkermansia_muciniphila, fill=infant_place_delivery)) + 
  labs (y="CLR transformed abundance", x="Place of delivery") + 
  ylim(-4, 13) +
  geom_sina(aes(color=infant_place_delivery), size=0.6) +
  geom_boxplot(aes(color=infant_place_delivery), outlier.shape = NA, alpha=0.5) + 
  scale_color_manual(name = "Place of delivery", 
                     labels=c('Home', 'Hospital', 'CS', 'VG'),
                     values=c("#FAE3D9", "#BBDED6", 'red', 'black')) +
  annotate("text", x=1.5, y=12, label="FDR=0.04,\nbeta hospital=-2.50817023") +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Place of delivery", 
                    labels=c('Home', 'Hospital'),
                    values=c("#FAE3D9", "#BBDED6")) 



vOTUs_phenos <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", RPKM_counts_VLP_filt_CLR, c("infant_place_delivery", 
                                                                                               "infant_ffq_feeding_mode_complex",
                                                                                               "infant_ever_never_breastfed"), 'time_as_covariate')

vOTUs_crossectional_solid_food <- linear_model_taxa(VLP_metadata[VLP_metadata$Timepoint=='M6',], "Short_sample_ID", RPKM_counts_VLP_filt_CLR, "infant_food_solid_intro_M6")

vOTUs_aggr_BacSp_phenos <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", vOTUs_assigned_BacSp_filt_CLR, c("infant_place_delivery", 
                                                                                                               "infant_ffq_feeding_mode_complex",
                                                                                                               "infant_ever_never_breastfed"), 'time_as_covariate')
#### FOR SUPPLEMENTARY TABLE ####
write.table(vOTUs_aggr_BacSp_phenos, '05.MANUSCRIPT/Supplementary_tables/MM_vOTUs_aggr_BacSp_RA_phenos_over_time_infant.txt', sep='\t', quote=F)


vOTUs_aggr_crossectional_solid_food <- linear_model_taxa(VLP_metadata[VLP_metadata$Timepoint=='M6',], "Short_sample_ID", vOTUs_assigned_BacSp_filt_CLR, "infant_food_solid_intro_M6")

place_of_delivery_phages <- vOTUs_assigned_BacSp_filt_CLR[,unique(vOTUs_aggr_BacSp_phenos[vOTUs_aggr_BacSp_phenos$P < 0.05 & vOTUs_aggr_BacSp_phenos$Pheno %in% c('infant_place_delivery'),]$Bug), drop=F]
colnames(place_of_delivery_phages) <- gsub('.*s__', '', colnames(place_of_delivery_phages))
place_of_delivery_phages <- merge(place_of_delivery_phages, VLP_metadata[,c("Timepoint", 'infant_place_delivery', 'infant_mode_delivery', 'Age_months', 'DNA_CONC', 'Clean_reads', 'NEXT_ID')], by='row.names')

pdf('./04.PLOTS/Phages_Akkermansia_muciniphila_infants_birthplace_no_CS.pdf', width=12/2.54, height=9/2.54)
ggplot(place_of_delivery_phages, aes(Timepoint, `Akkermansia muciniphila`, fill=infant_place_delivery)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Place of delivery", 
                    labels=c('Home', 'Hospital'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()

ggplot(place_of_delivery_phages, aes(infant_place_delivery, `Akkermansia muciniphila`, fill=infant_place_delivery)) + 
  labs (y="CLR transformed abundance", x="Place of delivery") + 
  ylim(-4, 13) +
  geom_sina(aes(color=infant_place_delivery), size=0.6) +
  geom_boxplot(aes(color=infant_place_delivery), outlier.shape = NA, alpha=0.5) + 
  scale_color_manual(name = "Place of delivery", 
                     labels=c('Home', 'Hospital', 'CS', 'VG'),
                     values=c("#FAE3D9", "#BBDED6", 'red', 'black')) +
  annotate("text", x=1.5, y=12, label="p-value=0.01,\nbeta hospital=-1.5") +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  scale_fill_manual(name = "Place of delivery", 
                    labels=c('Home', 'Hospital'),
                    values=c("#FAE3D9", "#BBDED6")) 

place_of_delivery_phages$Universal_fecal_ID <- gsub('V', '', place_of_delivery_phages$Row.names)
place_of_delivery_bacteria$Universal_fecal_ID <- gsub('B', '', place_of_delivery_bacteria$Row.names)
A_muciniphila_VH <- merge(place_of_delivery_bacteria[,c("Universal_fecal_ID", "Akkermansia_muciniphila", 
                                                        "Timepoint", "DNA_CONC", "NEXT_ID", "Clean_reads", 
                                                        "Age_months", "infant_place_delivery", "infant_mode_delivery")], place_of_delivery_phages[,c("Universal_fecal_ID", "Akkermansia muciniphila")], by="Universal_fecal_ID")
colnames(A_muciniphila_VH)[10] <- paste0('vOTUs_', colnames(A_muciniphila_VH)[10])
summary(lmer(`vOTUs_Akkermansia muciniphila` ~ Akkermansia_muciniphila + DNA_CONC + Clean_reads + infant_mode_delivery + infant_place_delivery + Age_months + (1|NEXT_ID) , REML = F,  data = A_muciniphila_VH))

ggplot(A_muciniphila_VH, aes(Akkermansia_muciniphila, `vOTUs_Akkermansia muciniphila`, color=infant_place_delivery)) + 
  geom_point() +
  geom_smooth(inherit.aes = F, aes(Akkermansia_muciniphila, `vOTUs_Akkermansia muciniphila`), method="lm") +
  annotate("text", x=4, y=7, label="p-value=5.1e-09\nbeta bacteria=2.6-01")+
  theme_bw()
  
host_assignment <- read.table('02.CLEAN_DATA/Host_prediction_to_genus_m90_refined_taxonomy_no_generalists.txt', sep='\t', header=T)

###### OUTPUT #####
write.table(species_phenos, '03a.RESULTS/TOP_PREVALENT_RA_SPECIES_phenos_mixed_models_all_results.txt', sep='\t', row.names=F, quote=F)
write.table(vOTUs_phenos, '03a.RESULTS/TOP_PREVALENT_RA_vOTUs_phenos_mixed_models_all_results.txt', sep='\t', row.names=F, quote=F)
write.table(vOTUs_aggr_BacSp_phenos, '03a.RESULTS/TOP_PREVALENT_RA_aggregated_vOTUs_at_species_phenos_mixed_models_all_results.txt', sep='\t', row.names=F, quote=F)

###### FOR VISUALIZATION #####
write.table(place_of_delivery_bacteria, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3G_A_muciniphila_place_of_birth.txt', sep='\t', row.names = F, quote=F)
write.table(place_of_delivery_phages, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3H_A_muciniphila_phages_place_of_birth.txt', sep='\t', row.names = F, quote=F)
write.table(A_muciniphila_VH, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3I_A_muciniphila_correlation_with_its_phages.txt', sep='\t', row.names = F, quote=F)
