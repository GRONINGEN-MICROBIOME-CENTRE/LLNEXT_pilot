setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we test the transmission enrichment in post vs pre
# birth maternal samples
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################

selected_viruses <- read.table("02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt", sep='\t', header=T)
selected_viruses <- selected_viruses[!is.na(selected_viruses$Transmitted_in_N_related),]
selected_viruses <- selected_viruses[selected_viruses$Transmitted_in_N_related!=0,]

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)

virus <- virus[names(virus) %in% selected_viruses$Virus]

selected_viruses$cutpoint_virus <- ifelse(selected_viruses$Youden_index >= selected_viruses$FDR_ipv_Youden, selected_viruses$FDR_ipv_Youden, selected_viruses$Youden_index)
selected_viruses$cutpoint_bacterium <- ifelse(selected_viruses$Host_Youden_index >= selected_viruses$Host_FDR_ipv_Youden, selected_viruses$Host_FDR_ipv_Youden, selected_viruses$Host_Youden_index)

virus <- lapply(virus, function(x) {
  x <- x/median( x[upper.tri(x)] )
})

#### exchanging distances values with 0 and 1 depending on the threshold for strain sharing event:
for (virusName in names(virus)) {
  
  virusN <- virus[[virusName]]
  
  virusN[virusN <= selected_viruses[selected_viruses$Virus==virusName,]$cutpoint_virus] <- 0
  
  virusN[virusN > selected_viruses[selected_viruses$Virus==virusName,]$cutpoint_virus] <- 1
  
  virus[[virusName]] <- virusN
}

metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/metadata_combined_for_exp.txt', sep='\t', header=T)
unique(metadata[metadata$Type=="Mother",]$Timepoint)
# testing:

preB_time <- c("P3", "P7", "B")
postB_time <- c("M1", "M2", "M3")


testing_enrichment_transmission <- list()


for ( virusName in names(virus) ) {
  
  virusN <- virus[[virusName]]
  virusN <- virusN[ grep('Infant', row.names(virusN)) , grep('Mother', colnames(virusN)) ]
  
  N_preB_transmission <- numeric()
  N_preB_copresence <- numeric()
  N_postB_transmission <- numeric()
  N_postB_copresence <- numeric()
  
  for (i in unique( gsub('_.*', '', colnames(virusN)) )  ) {
    
    if ( length( grep(i, row.names(virusN)) ) > 0 ) {
      
      A <- virusN[grep(i, row.names(virusN)), grep(i, colnames(virusN)), drop=F]
      
      N_preB_transmission[i] <- sum(unlist(A[, grep(paste(preB_time,collapse="|"), 
                                                    colnames(A), value=TRUE), drop=F ])==0)
      N_preB_copresence[i] <- length(unlist(A[, grep(paste(preB_time,collapse="|"), 
                                                     colnames(A), value=TRUE), drop=F ]))
      
      N_postB_transmission[i] <- sum(unlist(A[, grep(paste(postB_time,collapse="|"), 
                                                     colnames(A), value=TRUE), drop=F ])==0)
      N_postB_copresence[i] <- length(unlist(A[, grep(paste(postB_time,collapse="|"), 
                                                      colnames(A), value=TRUE), drop=F ]))
    }
    
    
  }
  
  testing_enrichment_transmission[[virusName]] <- matrix( c(sum(N_preB_transmission), # transmitted from preB
                                                            sum(N_postB_transmission), # transmitted from postB
                                                            sum(N_preB_copresence) -  sum(N_preB_transmission), # not transmitted from preB
                                                            sum(N_postB_copresence) -  sum(N_postB_transmission) ), # not transmitted from postB
                                                          nrow=2,
                                                          dimnames = list(c('preB', 'postB'),
                                                                          c('Transmitted', 'Not-transmitted')
                                                                          )) 
  
}

colnames(virusN)



selected_viruses <- selected_viruses[!duplicated(selected_viruses$Virus),]
selected_viruses <- selected_viruses[selected_viruses$Transmission_enriched_in_related=="YES",]
selected_viruses$Different_between_preB_postB <- NA

for (virusName in selected_viruses$Virus) {
  
  fisher_test <- fisher.test(testing_enrichment_transmission[[virusName]])
  
  selected_viruses[selected_viruses$Virus==virusName, "Different_between_preB_postB"] <- fisher_test$p.value
  
}

selected_viruses$FDR_difference_freq_time <- p.adjust(selected_viruses$Different_between_preB_postB)
View(selected_viruses[selected_viruses$FDR_difference_freq_time< 0.05,])

fisher_test <- fisher.test(testing_enrichment_transmission[["LN_7B08_VL_396_NODE_588_length_5274_cov_19.581337"]])

mosaicplot(testing_enrichment_transmission[["LN_7B08_VL_396_NODE_588_length_5274_cov_19.581337"]], color=T)

keeper <- selected_viruses
#### FOR SUPPLEMENTARY TABLE ####
write.table(keeper, '05.MANUSCRIPT/Supplementary_tables/Virus_sharing_testing_enrichment_pre_postbirth.txt', sep='\t', quote=F, row.names=F)


MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_', MGS_metadata$source, '_', MGS_metadata$Timepoint)

selected_viruses <- read.table("02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt", sep='\t', header=T)
selected_viruses <- selected_viruses[!is.na(selected_viruses$Transmitted_in_N_related),]
selected_viruses <- selected_viruses[selected_viruses$Transmitted_in_N_related!=0,]

selected_viruses <- selected_viruses[!is.na(selected_viruses$Host_Transmission_enriched_in_related),]
selected_viruses <- selected_viruses[selected_viruses$Host_Transmission_enriched_in_related=="YES",]
selected_viruses <- selected_viruses[!duplicated(selected_viruses$Host_SGB),]


folder = "01.RAW_DATA/strainphlan_4_distance_matrix/"  
file_list = list.files(path=folder, pattern="*.csv")  
bacterium=lapply(paste0(folder, file_list), function(x) read.csv(x, header=T))
names(bacterium) <- gsub("_dmat_.csv", '', file_list)
bacterium <- bacterium[names(bacterium) %in% selected_viruses$Host_SGB]
bacterium <- lapply(bacterium, function(x) {
  x <- x/median( x[upper.tri(x)], na.rm=T )
})

## since my parsing scripts have certain ID format: 
bacterium <- lapply(bacterium, function(x) {
  row.names(x) <- substr(row.names(x), 0, 12) # first remove data processing artifacts
  colnames(x) <- substr(colnames(x), 0, 12) 
  
  row.names(x) <- MGS_metadata$Alex_ID[match(row.names(x), MGS_metadata$NG_ID)]
  colnames(x) <- MGS_metadata$Alex_ID[match(colnames(x), MGS_metadata$NG_ID)] # then change the sequencing IDs to IDs for parsing
  x
})

#### exchanging distances values with 0 and 1 depending on the threshold for strain sharing event:
for (bacteriumName in names(bacterium)) {
  
  bacteriumN <- bacterium[[bacteriumName]]
  
  bacteriumN[bacteriumN <= selected_viruses[selected_viruses$Host_SGB==bacteriumName,]$cutpoint_bacterium] <- 0
  
  bacteriumN[bacteriumN > selected_viruses[selected_viruses$Host_SGB==bacteriumName,]$cutpoint_bacterium] <- 1
  
  bacterium[[bacteriumName]] <- bacteriumN
}

testing_enrichment_transmission_bacterium <- list()


selected_viruses$Different_between_preB_postB_bact <- NA

for ( bacteriumName in names(bacterium) ) {
  
  bacteriumN <- bacterium[[bacteriumName]]
  bacteriumN <- bacteriumN[ grep('Infant', row.names(bacteriumN)) , grep('Mother', colnames(bacteriumN)) ]
  
  N_preB_transmission <- numeric()
  N_preB_copresence <- numeric()
  N_postB_transmission <- numeric()
  N_postB_copresence <- numeric()
  
  for (i in unique( gsub('_.*', '', colnames(bacteriumN)) )  ) {
    
    if ( length( grep(i, row.names(bacteriumN)) ) > 0 ) {
      
      A <- bacteriumN[grep(i, row.names(bacteriumN)), grep(i, colnames(bacteriumN)), drop=F]
      
      N_preB_transmission[i] <- sum(unlist(A[, grep(paste(preB_time,collapse="|"), 
                                                    colnames(A), value=TRUE), drop=F ])==0, na.rm = T)
      N_preB_copresence[i] <- length(unlist(A[, grep(paste(preB_time,collapse="|"), 
                                                     colnames(A), value=TRUE), drop=F ]))
      
      N_postB_transmission[i] <- sum(unlist(A[, grep(paste(postB_time,collapse="|"), 
                                                     colnames(A), value=TRUE), drop=F ])==0, na.rm = T)
      N_postB_copresence[i] <- length(unlist(A[, grep(paste(postB_time,collapse="|"), 
                                                      colnames(A), value=TRUE), drop=F ]))
    }
    
    
  }
  
  testing_enrichment_transmission_bacterium[[bacteriumName]] <- matrix( c(sum(N_preB_transmission), # transmitted from preB
                                                            sum(N_postB_transmission), # transmitted from postB
                                                            sum(N_preB_copresence) -  sum(N_preB_transmission), # not transmitted from preB
                                                            sum(N_postB_copresence) -  sum(N_postB_transmission) ), # not transmitted from postB
                                                          nrow=2,
                                                          dimnames = list(c('preB', 'postB'),
                                                                          c('Transmitted', 'Not-transmitted')
                                                          )) 
  
}

for (bacteriumName in selected_viruses$Host_SGB) {
  
  fisher_test <- fisher.test(testing_enrichment_transmission_bacterium[[bacteriumName]])
  
  selected_viruses[selected_viruses$Host_SGB==bacteriumName, "Different_between_preB_postB_bact"] <- fisher_test$p.value
  
}

mosaicplot(testing_enrichment_transmission_bacterium[["SGB1934"]], color=T)
selected_viruses$FDR_Different_between_preB_postB_bact <- p.adjust(selected_viruses$Different_between_preB_postB_bact, method = "BH")
write.table(selected_viruses, '05.MANUSCRIPT/Supplementary_tables/BacStr_sharing_testing_enrichment_pre_postbirth.txt', sep='\t', quote=F, row.names=F)

