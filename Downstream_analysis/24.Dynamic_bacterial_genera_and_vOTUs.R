setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore changes in abundance of gut bacterial
# genera and aggregated vOTUs at the level of bacterial genus
# in the developing infant gut 
#############################################################
# Authors: Trishla Sinha & Sanzhima Garmaeva

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
library(dplyr)
library(tibble)
library(lme4)
library(RLRsim)
library(lmerTest)

library(reshape2)
library(ggplot2)
library(ggforce)

library(tidyverse)
library(ape)
library(ggtree)
##############################
# Input data
##############################
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata <- MGS_metadata[MGS_metadata$Type=='Infant',]
MGS_metadata$Timepoint <- factor(MGS_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M9", "M12"), ordered = T)
row.names(MGS_metadata) <- MGS_metadata$Short_sample_ID_bact

genera <- read.table('02.CLEAN_DATA/Microbiome_genera_unfiltred.txt', sep='\t', header=T)
genera <- genera[,MGS_metadata$Short_sample_ID_bact]
# filtering for presence in more than 10% of microbiome samples:
genera_filt <- genera[(rowSums(genera!=0) > 0.10*ncol(genera)), ]
genera_filt <- as.data.frame(t(genera_filt))
# CLR-transformation
my_pseudocount_normal=min(genera_filt[genera_filt!=0])/2
genera_filt_CLR<-decostand(genera_filt, "clr", pseudocount=my_pseudocount_normal)

microbiome_genera <- gsub(".*g__", "", row.names(genera))

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata <- VLP_metadata[VLP_metadata$Type=='Infant',]
VLP_metadata$Timepoint <- factor(VLP_metadata$Timepoint, levels = c("M1", "M2", "M3",
                                                                    "M6", "M12"), ordered = T)
row.names(VLP_metadata) <- VLP_metadata$Short_sample_ID

vOTUS_assigned_BacGen <- read.table('02.CLEAN_DATA/RPKM_counts_aggregated_bacterial_host_genus.txt', sep='\t', header=T)
row.names(vOTUS_assigned_BacGen) <- vOTUS_assigned_BacGen$Host.genus
vOTUS_assigned_BacGen <- vOTUS_assigned_BacGen[,VLP_metadata$Short_sample_ID]
vOTUS_assigned_BacGen_filt <- vOTUS_assigned_BacGen[(rowSums(vOTUS_assigned_BacGen!=0) > 0.10*ncol(vOTUS_assigned_BacGen)),]
vOTUS_assigned_BacGen_filt <- as.data.frame(t(vOTUS_assigned_BacGen_filt))
#CLR-transformation
my_pseudocount_normal=min(vOTUS_assigned_BacGen_filt[vOTUS_assigned_BacGen_filt!=0])/2
vOTUS_assigned_BacGen_filt_CLR <- decostand(vOTUS_assigned_BacGen_filt, "clr", pseudocount=my_pseudocount_normal)

colors_Fig2g <- read.table("02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig2GH_gen_palette_colors.txt", sep='\t', header=T) 
##############################
# ANALYSIS
##############################
# Choose only those vOTU-host-genus aggregates that are present in prevalent bacterial genera:
vOTUsgen_correlate <- vOTUS_assigned_BacGen_filt_CLR[,gsub('.*g__','',colnames(vOTUS_assigned_BacGen_filt_CLR)) %in% gsub('.*g__','',colnames(genera_filt_CLR))]
row.names(vOTUsgen_correlate) <- gsub('V', '', row.names(vOTUsgen_correlate) )
# Choose only those bacterial genera that are present among vOTU-host-genus aggregates:
bacgen_correlate <- genera_filt_CLR[,gsub('.*g__','',colnames(genera_filt_CLR)) %in% gsub('.*g__','',colnames(vOTUS_assigned_BacGen_filt_CLR))]
bacgen_correlate$Universal_fecal_ID <- gsub('B', '', row.names(bacgen_correlate) )
bac_meta <- merge(bacgen_correlate, MGS_metadata[,c("Universal_fecal_ID", "DNA_CONC", "Clean_reads", 
                                                    "infant_mode_delivery", "Age_months", "Individual_ID")],  
                  by='Universal_fecal_ID')

virus_host_interaction <- mixed_models_taxa(bac_meta, "Universal_fecal_ID", vOTUsgen_correlate, grep('k__', colnames(bacgen_correlate), value=T), 'time_as_covariate')
virus_host_interaction$Feature <- gsub('`','', gsub('.*g__', '', virus_host_interaction$Feature))
virus_host_interaction$Bug <- paste0('vOTUs_', gsub('.*g__', '', virus_host_interaction$Bug))

# actual testing
#virus_host_interaction <- mixed_models_VH_interaction(VLP_metadata, "Universal_fecal_ID", host_virus)
virus_host_interaction$Significane_level <- NA
virus_host_interaction[virus_host_interaction$FDR > 0.05,]$Significane_level <- ""
virus_host_interaction[virus_host_interaction$FDR <= 0.05,]$Significane_level <- "*"
virus_host_interaction[virus_host_interaction$FDR <= 0.01,]$Significane_level <- "**"
virus_host_interaction[virus_host_interaction$FDR <= 0.001,]$Significane_level <- "***"

# filtering for those prevalent in at least 20% samples
virus_host_interaction_0.2perc <- virus_host_interaction[virus_host_interaction$Feature %in% gsub('.*g__', '', row.names(vOTUS_assigned_BacGen[(rowSums(vOTUS_assigned_BacGen!=0) > 0.20*ncol(vOTUS_assigned_BacGen)),])),]
virus_host_interaction_0.2perc <- virus_host_interaction_0.2perc[order(virus_host_interaction_0.2perc$Pheno),]
virus_host_interaction_0.2perc$Feature <- factor(virus_host_interaction_0.2perc$Feature, levels=unique(virus_host_interaction_0.2perc$Feature), ordered = T)
virus_host_interaction_0.2perc <- virus_host_interaction_0.2perc[virus_host_interaction_0.2perc$Bug %in% paste0('vOTUs_', gsub('.*g__', '', row.names(vOTUS_assigned_BacGen[(rowSums(vOTUS_assigned_BacGen!=0) > 0.20*ncol(vOTUS_assigned_BacGen)),]))),]
virus_host_interaction_0.2perc$diag <- virus_host_interaction_0.2perc$Feature == gsub('vOTUs_', '', virus_host_interaction_0.2perc$Bug)
virus_host_interaction_0.2perc$Bug <- factor(virus_host_interaction_0.2perc$Bug, levels=paste0('vOTUs_',unique(virus_host_interaction_0.2perc$Feature)), ordered = T)

pdf('./04.PLOTS/Prevalent_at_0.2_vOTUs_aggregates_and_their_hosts.pdf', width=21/2.54, height=21/2.54)
ggplot(virus_host_interaction_0.2perc, aes(Feature, Bug, fill = Estimate))+
  labs(x='CLR-t abundance of bacterial genus', y='CLR-t abundance of vOTUs-host-genus aggregate') +
  geom_tile(aes(color=diag), width = 0.98, height = 0.98) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Effect size") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text( vjust = 1, 
                                   size = 10, hjust = 1)) +
  geom_text(aes(label=Significane_level), color="black", size=3) + 
  scale_color_manual(guide = FALSE, values = c("TRUE" = "black", "FALSE" = "white")) +
  coord_flip()
dev.off()

dat <- as.data.frame(unique(virus_host_interaction_0.2perc$Pheno))
dat <- separate(dat, sep = '\\|', col = "unique(virus_host_interaction_0.2perc$Pheno)", into = c('Kingdom', 'Phylum', 'Class', 
                                                                                        'Order', 'Family', 'Genus'))
dat <- lapply( dat, function(x){
  gsub('.*.__', '', x)
}  )

dat <- data.frame(lapply(dat[sapply(dat, is.character)], 
                            as.factor))


taxa <- as.phylo(~Kingdom/Phylum/Class/Order/Family/Genus, data = dat)

pdf('./04.PLOTS/Phylo.pdf', width=5/2.54, height=21/2.54)
plot(taxa, type = "cladogram", align.tip.label = T)
dev.off()

# Dynamic bacterial genera (association of bacterial genus abundance with numeric time):

# association with timepoints:
genera_dynamics <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", genera_filt_CLR, "Age_months", "don't consider time as a covariate")
genera_dynamics$Genus <- gsub(".*g__", "", genera_dynamics$Bug)
genera_dynamics <- genera_dynamics[genera_dynamics$Genus %in% unique(virus_host_interaction_0.2perc$Feature), ]
genera_dynamics <- genera_dynamics[genera_dynamics$FDR <= 0.05,]
# association with phenotypes corrected for time
genera_phenos_UPD <- mixed_models_taxa(MGS_metadata, "Short_sample_ID_bact", genera_filt_CLR, c("infant_place_delivery", 
                                                                                                "infant_ffq_feeding_mode_complex"), 'time_as_covariate')

# Dynamic vOTUs aggregates at the host genus (association of vOTUs aggregates abundance with numeric time):
vOTUs_aggr_gen_dynamics <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", vOTUS_assigned_BacGen_filt_CLR, "Age_months", "don't consider time as a covariate")
vOTUs_aggr_gen_dynamics$Genus <- gsub(".*g__", "", vOTUs_aggr_gen_dynamics$Bug)
##### FOR SUPPLEMENTARY TABLE ####
write.table(vOTUs_aggr_gen_dynamics, "05.MANUSCRIPT/Supplementary_tables/DYNAMIC_TOP_PREVALENT_RA_vOTUs_aggregates_FIRST_YEAR_NEXT_PILOT_07_07_2023.txt", sep='\t', quote=F)

vOTUs_aggr_gen_dynamics <- vOTUs_aggr_gen_dynamics[vOTUs_aggr_gen_dynamics$Genus %in% unique(virus_host_interaction_0.2perc$Feature), ]
vOTUs_aggr_gen_dynamics <- vOTUs_aggr_gen_dynamics[vOTUs_aggr_gen_dynamics$FDR <= 0.05,]


# association vOTUs aggregates at the host genus level with phenotypes corrected for time
vOTUs_gen_phenos <- mixed_models_taxa(VLP_metadata, "Short_sample_ID", vOTUS_assigned_BacGen_filt_CLR, c("infant_place_delivery", 
                                                                                                         "infant_ffq_feeding_mode_complex"), 'time_as_covariate')

dynamic_to_visualize <- virus_host_interaction_0.2perc[virus_host_interaction_0.2perc$diag==T,]
dynamic_to_visualize <- dynamic_to_visualize[order(dynamic_to_visualize$FDR),]
dynamic_to_visualize <- dynamic_to_visualize[dynamic_to_visualize$Feature %in% genera_dynamics$Genus,]
dynamic_to_visualize <- dynamic_to_visualize[dynamic_to_visualize$Feature %in% vOTUs_aggr_gen_dynamics$Genus,]
dynamic_to_visualize <- dynamic_to_visualize[dynamic_to_visualize$FDR <= 0.05,]

# plotting top associates prevalent in 20% and associated with time:

# 1. vOTUs: 
dynamic_vOTUs_plot <- vOTUS_assigned_BacGen_filt_CLR[,gsub('.*g__', '', colnames(vOTUS_assigned_BacGen_filt_CLR)) %in% dynamic_to_visualize$Feature]
dynamic_vOTUs_plot <- merge(dynamic_vOTUs_plot, VLP_metadata[,c("Timepoint", "NEXT_ID")], by='row.names')
dynamic_vOTUs_plot_melt <- melt(dynamic_vOTUs_plot)

dynamic_vOTUs_plot_melt <- dynamic_vOTUs_plot_melt[order(dynamic_vOTUs_plot_melt$NEXT_ID, dynamic_vOTUs_plot_melt$Timepoint),]
dynamic_vOTUs_plot_melt <- dynamic_vOTUs_plot_melt %>% 
  arrange(NEXT_ID, Timepoint)

dynamic_vOTUs_plot_melt$variable <- gsub('.*g__','',dynamic_vOTUs_plot_melt$variable)

trial_VLP <-dynamic_vOTUs_plot_melt  %>% group_by(variable, Timepoint) %>% 
  summarise(mean = mean(value), n = n())

# COLOR SCHEME:
colors_Fig2g <- colors_Fig2g[colors_Fig2g$Genus %in% unique(dynamic_vOTUs_plot_melt$variable),]

colors <- data.frame(sort(unique(dynamic_vOTUs_plot_melt$variable)),
                     c("#2f357c","#FF0000", "#ff0080", 
                       "#FFA500", "#0000FF", "#b0799a",
                       "#800000", "#FFC0CB", "#008000", 
                       "#008080", "#DC143C", "#bf3729",
                       "#40E0D0",  "#898e9f"))

colnames(colors) <- c("Genus", "Color")
show_col(colors$Color)


pdf('./04.PLOTS/Dynamic_phages_aggregated_by_bacterial_genera_infants.pdf', width=16/2.54, height=14/2.54)
ggplot(trial_VLP, aes(x = Timepoint, y = mean, color = variable, group = variable)) +
  geom_point(size=2) +
  geom_line(size=1.5, alpha=0.8) +
  scale_color_manual(values = colors$Color) +
  ggtitle("")+
  theme_bw()+
  labs(x="Timepoint", y = "Mean CLR Tranformed abundance", color='Bacterial\nhost genera')+
  theme(
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12,  hjust = 1), 
    legend.title = element_text(color="black", size=12, face="bold"))
dev.off()

# 2. Bacterial genera: 

dynamic_plot <- genera_filt_CLR[,gsub('.*g__', '', colnames(genera_filt_CLR)) %in% dynamic_to_visualize$Feature]
dynamic_plot <- merge(dynamic_plot, MGS_metadata[,c("Timepoint", "NEXT_ID")], by='row.names')
dynamic_plot_melt <- melt(dynamic_plot)

dynamic_plot_melt <- dynamic_plot_melt[order(dynamic_plot_melt$NEXT_ID, dynamic_plot_melt$Timepoint),]
dynamic_plot_melt <- dynamic_plot_melt %>% 
  arrange(NEXT_ID, Timepoint)

dynamic_plot_melt$variable <- gsub('.*g__','',dynamic_plot_melt$variable)
  
trial <- dynamic_plot_melt  %>% group_by(variable, Timepoint) %>% 
  summarise(mean = mean(value), n = n())

pdf('./04.PLOTS/Dynamic_bacterial_genera_infants_24_05_2023.pdf', width=16/2.54, height=14/2.54)
ggplot(trial, aes(x = Timepoint, y = mean, color = variable, group = variable)) +
  geom_point(size=2) +
  geom_line(size=1.5, alpha=0.8) +
  scale_color_manual(values = colors$Color) +
  ggtitle("")+
  theme_bw()+
  labs(x="Timepoint", y = "Mean CLR Tranformed abundance", color='Bacterial Genera')+
  theme(
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    axis.text.y = element_text(face="bold", size=10),
    axis.text.x = element_text(size=12, hjust = 1), 
    legend.title = element_text(color="black", size=12, face="bold"))
dev.off()


# correlations:
vOTUsgen_correlate$Universal_fecal_ID <- row.names(vOTUsgen_correlate)
correlations_to_plot <- merge(bacgen_correlate, vOTUsgen_correlate, by='Universal_fecal_ID')
row.names(correlations_to_plot) <- correlations_to_plot$Universal_fecal_ID
# filtering only those that have a FDR-significant correlation:
FDR_significant_correlation <-  virus_host_interaction_0.2perc[virus_host_interaction_0.2perc$diag==T & virus_host_interaction_0.2perc$FDR <= 0.05,"Feature"]

correlations_to_plot <- correlations_to_plot[, gsub('.*g__', '', colnames(correlations_to_plot)) %in% FDR_significant_correlation]
correlations_to_plot$Timepoint <- VLP_metadata$Timepoint[match(row.names(correlations_to_plot), VLP_metadata$Universal_fecal_ID)]
colnames(correlations_to_plot)[grep('k__', colnames(correlations_to_plot))] <- gsub('.*g__','',colnames(correlations_to_plot)[grep('k__', colnames(correlations_to_plot))])
colnames(correlations_to_plot)[grep('d__', colnames(correlations_to_plot))] <- paste0('vOTUs_', gsub('.*g__','',colnames(correlations_to_plot)[grep('d__', colnames(correlations_to_plot))]))

ggplot(correlations_to_plot, aes(Klebsiella, vOTUs_Klebsiella, color=Timepoint)) +
  geom_point() + 
  labs(x="CLR-transformed abundance of bacterial genus", y="CLR-transformed abundance of vOTUs-aggregate") +
  ggtitle("Klebsiella and their predicted phages\nover time") +
  annotate(geom = "text", x = 3, y=6, label="FDR=2.2e-06\nEffect size=0.3") +
  geom_smooth(inherit.aes = F, aes(Klebsiella, vOTUs_Klebsiella), method = "lm") +
  theme_bw() + 
  theme(plot.title = element_text(face="bold", hjust=0.5))

# selecting genera dependant on phenotypes:
genera_phenos_FDR <- genera_phenos_UPD[genera_phenos_UPD$FDR < 0.05,]

# selected genera dependent on feeding mode:
feeding_plot <- genera_filt_CLR[,unique(genera_phenos_FDR[genera_phenos_FDR$Pheno=='infant_ffq_feeding_mode_complex',]$Bug)]
colnames(feeding_plot) <- gsub('.*g__', '', colnames(feeding_plot))
feeding_plot <- merge(feeding_plot, MGS_metadata[,c("Timepoint", 'infant_ffq_feeding_mode_complex')], by='row.names')
feeding_plot <- feeding_plot[!is.na(feeding_plot$infant_ffq_feeding_mode_complex),]

pdf('./04.PLOTS/Genus_Cutibacterium_infants_feeding.pdf', width=13/2.54, height=9/2.54)
ggplot(feeding_plot, aes(Timepoint, Cutibacterium, fill=infant_ffq_feeding_mode_complex)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  ggtitle("g__Cutibacterium") +
  scale_fill_manual(name = "Feeding mode", 
                    labels=c('BF', 'FF', 'MF'),
                    values=c("#143F6B", "#F55353", "#FEB139"))
dev.off()

# selected genus dependent on place of delivery:

birth_plot <- genera_filt_CLR[,unique(genera_phenos_FDR[genera_phenos_FDR$Pheno %in% c('infant_place_delivery'),]$Bug), drop=F]
colnames(birth_plot) <- gsub('.*g__', '', colnames(birth_plot))
birth_plot <- merge(birth_plot, MGS_metadata[,c("Timepoint", 'infant_place_delivery', 'infant_mode_delivery', 'Age_months', 'DNA_CONC', 'Clean_reads', 'NEXT_ID')], by='row.names')

# for additional analysis for presence of Akkermansia
Akkermansia_binary <- genera[grep('Akkermansia', row.names(genera)),] 
Akkermansia_binary[Akkermansia_binary!=0] <- 1
Akkermansia_binary <- as.data.frame(t(Akkermansia_binary))
colnames(Akkermansia_binary) <- 'Akkermansia_binary'
Akkermansia_binary$Row.names <- row.names(Akkermansia_binary)

birth_plot <- merge(birth_plot, Akkermansia_binary, by='Row.names')

place_mod0 <- lmer(Akkermansia ~ Age_months + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
place_mod1  = lmer(Akkermansia ~ Age_months + DNA_CONC + Clean_reads + infant_place_delivery + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
BIC(place_mod0, place_mod1)
summary(place_mod1) # -2.528e+00  7.217e-01  2.787e+01  -3.502  0.00157 **

place_mod2 <- lmer(Akkermansia_binary ~ Age_months + DNA_CONC + Clean_reads + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
place_mod3  = lmer(Akkermansia_binary ~ Age_months + DNA_CONC + Clean_reads + infant_place_delivery + (1|NEXT_ID), REML = F, data = birth_plot[birth_plot$infant_mode_delivery=='VG',])
BIC(place_mod2, place_mod3)
summary(place_mod3) # infant_place_deliveryhospital -2.424e-01  7.383e-02  2.806e+01  -3.283  0.00275 **

pdf('./04.PLOTS/Genus_Akkermansia_infants_birthplace_no_CS.pdf', width=12/2.54, height=9/2.54)
ggplot(birth_plot[birth_plot$infant_mode_delivery=='VG',], aes(Timepoint, Akkermansia, fill=infant_place_delivery)) + 
  labs (y="CLR transformed abundance", x="Timepoint") + 
  geom_boxplot(outlier.shape = NA, alpha=0.5) + 
  geom_sina(colour='#303841', size=0.6) +
  theme_bw()+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=16,face="bold"),
        strip.text.x = element_text(size = 12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, face="bold")) +
  ggtitle("g__Akkermansia") +
  scale_fill_manual(name = "Place of delivery", 
                    labels=c('Home', 'Hospital'),
                    values=c("#FAE3D9", "#BBDED6"))
dev.off()

###### OUTPUT #####
write.table(dynamic_genera, '03a.RESULTS/DYNAMIC_TOP_PREVALENT_RA_GENERA_FIRST_YEAR_NEXT_PILOT_07_07_2023.txt', sep='\t', row.names = F, quote=F)
write.table(genera_phenos_UPD, '03a.RESULTS/TOP_PREVALENT_RA_GENERA_phenos_mixed_models_all_results_07_07_2023.txt', sep='\t', row.names = F, quote=F)
write.table(vOTUs_gen_phenos, '03a.RESULTS/TOP_PREVALENT_vOTUs_aggr_Bac_Gen_phenos_mixed_models_all_results_07_07_2023.txt', sep='\t', row.names = F, quote=F)

###### FOR VISUALIZATION #####
write.table(virus_host_interaction_0.2perc, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3A_mixed_model_outcome_vOTUs-host-genus_aggregates_at_0.2_prevalence_vs_bacgen.txt', sep='\t', row.names=F, quote=F)
write.table(correlations_to_plot, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3B_FDR_correlated_bacteria_vs_its_predicted_phages.txt', sep='\t', quote=F)
write.table(trial_VLP, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3A_mean_vOTUs_aggregate_genus_p_time_infant.txt', sep='\t', row.names=F, quote=F)
write.table(colors, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3AB_color_scheme.txt', sep='\t', row.names=F)
write.table(trial, '02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/SFig3B_mean_bacterial_genus_p_time_infant.txt', sep='\t', row.names=F, quote=F)
