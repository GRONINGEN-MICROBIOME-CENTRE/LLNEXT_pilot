setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

library(vegan)
library(lme4)
library(RLRsim)
library(lmerTest)
VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T)

MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T)
VLP_metadata <- VLP_metadata[VLP_metadata$Universal_fecal_ID %in% MGS_metadata$Universal_fecal_ID,]
MGS_metadata <- MGS_metadata[MGS_metadata$Universal_fecal_ID %in% VLP_metadata$Universal_fecal_ID,]

VLP_metadata <- VLP_metadata[order(VLP_metadata$Universal_fecal_ID),]
MGS_metadata <- MGS_metadata[order(MGS_metadata$Universal_fecal_ID),]


RPKM_counts_VLP <- read.table("02.CLEAN_DATA/RPKM_counts_VLP.txt", sep='\t', header=T)
RPKM_counts_VLP <- RPKM_counts_VLP[, colnames(RPKM_counts_VLP) %in% VLP_metadata[VLP_metadata$Universal_fecal_ID %in% MGS_metadata$Universal_fecal_ID,]$Short_sample_ID]
RPKM_counts_VLP <- RPKM_counts_VLP[,VLP_metadata$Short_sample_ID]

RPKM_counts_MGS <- read.table("02.CLEAN_DATA/RPKM_counts_MGS.txt", sep='\t', header=T)
RPKM_counts_MGS <- RPKM_counts_MGS[, colnames(RPKM_counts_MGS) %in% MGS_metadata[MGS_metadata$Universal_fecal_ID %in% VLP_metadata$Universal_fecal_ID,]$Short_sample_ID]
RPKM_counts_MGS <- RPKM_counts_MGS[,MGS_metadata$Short_sample_ID]

RPKM_counts_combo <- merge(RPKM_counts_VLP, RPKM_counts_MGS, by='row.names', all = T)
row.names(RPKM_counts_combo) <- RPKM_counts_combo$Row.names
RPKM_counts_combo$Row.names <- NULL
RPKM_counts_combo[is.na(RPKM_counts_combo)] <- 0
RPKM_counts_combo <- RPKM_counts_combo[rowSums(RPKM_counts_combo)>0,]

dist = vegdist(as.data.frame(t(RPKM_counts_combo)), method = "bray")
dist_matrix <- as.matrix(dist)

see <- dist_matrix[1:204, 205:408]
median(diag(see))
round(sd(diag(see)), 2)

unsee <- see
for (i in gsub('FAM','',unique(VLP_metadata$FAM_ID))) {
  
  unsee[substr(colnames(unsee), 4, 7)==i, substr(row.names(unsee), 4, 7)==i] <- NA
  
}
round(median(unname(unlist(unsee)), na.rm=T), 2)
round(sd(unname(unlist(unsee)), na.rm=T), 2)

summary(lm(MGS_metadata$viral_alpha_diversity_MGS ~ VLP_metadata$N_timepoints))

type_mod0 <- lm(MGS_metadata$viral_richness_MGS ~ VLP_metadata$N_timepoints + MGS_metadata$DNA_CONC + MGS_metadata$Clean_reads)
type_mod1  = lmer(MGS_metadata$viral_richness_MGS ~ VLP_metadata$N_timepoints + MGS_metadata$DNA_CONC + MGS_metadata$Clean_reads + (1|MGS_metadata$NEXT_ID), REML = F, data = VLP_metadata)
BIC(type_mod0, type_mod1)
exactLRT(type_mod1,type_mod0)
summary(type_mod1) #Type 2.525e+00  3.517e-01 6.182e+01   7.178 1.07e-09 ***

type_mod0 <- lm(VLP_metadata$viral_alpha_diversity ~ VLP_metadata$N_timepoints + VLP_metadata$DNA_CONC + VLP_metadata$Clean_reads)
type_mod1  = lmer(VLP_metadata$viral_alpha_diversity ~ VLP_metadata$N_timepoints + VLP_metadata$DNA_CONC + VLP_metadata$Clean_reads + (1|VLP_metadata$NEXT_ID), REML = F, data = VLP_metadata)

real <- wilcox.test(diag(see), unname(unlist(unsee)))
real$p.value
p_value_perm <- list()

# loop over all bacteria

  
  # creating a vector for storing F-statistics for the bacterium[[n]]
  p_value <- as.numeric(c())
  
  # loop over 1000 permutations
  for (i in 1:1000) {
    
    # randomly permuting the samples:
    # first get the random permutation of columns
    sample_permut = sample( 1:ncol(see) )
    # second get the random permutation of columns
    sample_permut_row = sample( 1:nrow(see) )
    
    # reorder the table of distances for the chosen bacterium according to the permutation
    #perm_table = bacteriumN[sample_permut,sample_permut]
    perm_table = see[sample_permut_row,sample_permut]
    # rename columns
    colnames(perm_table) = colnames(see)
    row.names(perm_table) = row.names(see)
  
    unsee_perm_table <- perm_table
    for (k in gsub('FAM','',unique(VLP_metadata$FAM_ID))) {
      
      unsee_perm_table[substr(colnames(unsee_perm_table), 4, 7)==k, substr(row.names(unsee_perm_table), 4, 7)==k] <- NA
      
    }
    
    wilcox_perm=wilcox.test(diag(perm_table), unname(unlist(unsee_perm_table)), alternative='less', paired=F)
    
      # storing p-value for this permutation
      p_value[i] <- wilcox_perm$p.value
  }      

sum(p_value > real$p.value)

