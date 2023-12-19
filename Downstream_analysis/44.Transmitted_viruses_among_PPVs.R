setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')


#############################################################
# Here we explore if transmitted viruses are the members of 
# infant PPVs
#############################################################

##############################
# Functions
##############################
source("03.SCRIPTS/NEXT_pilot_FUP_downstream/stability_functions.R")
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
metadata$Alex_ID <- paste0(metadata$FAM_ID, '_', metadata$Type, '_', substr(metadata$Short_sample_ID, 1,1), '_', metadata$source, '_', metadata$Timepoint)

VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header = T)
VLP_metadata$Alex_ID_short <- paste0(VLP_metadata$FAM_ID, '_', VLP_metadata$Type, '_', substr(VLP_metadata$Short_sample_ID, 1,1) )
VLP_metadata <- VLP_metadata[VLP_metadata$Type=='Infant',]
VLP_metadata <- VLP_metadata[VLP_metadata$Individual_ID %in% names(PPVs),]

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)

vOTUs_infants <- RPKM_counts_VLP[,VLP_metadata[VLP_metadata$Type=='Infant',]$Short_sample_ID]
vOTUs_infants <- vOTUs_infants[rowSums(vOTUs_infants)!=0,]
##############################
# ANALYSIS
##############################
transmission_p_infant <- list()

for (virusName in names(virus)) {
 
  virusN <- virus[[virusName]]
  
  virusN <- virusN[grep('Mother', row.names(virusN)), grep('Infant', colnames(virusN))]
  
  for ( infant in unique(substr(metadata[metadata$Type=='Infant',]$Alex_ID, 1, 16)) ) {
    
    FAM <- substr(infant, 1, 7)
    
    if ( length(grep(FAM, row.names(virusN)))!=0 & length(grep(infant, colnames(virusN)))!=0 ) {
      
      infant_strains <- virusN[grep(FAM, row.names(virusN)),grep(infant, colnames(virusN)), drop=F]
      
      if ( any(infant_strains==0) ) {
        
        transmission_p_infant[[infant]] <- c(transmission_p_infant[[infant]], virusName)
        
      }
      
    }
    
  }
  
}


p_vir_frac_infants <- personal_biome_fraction(VLP_metadata[VLP_metadata$Type=='Infant' & VLP_metadata$N_timepoints>2,], vOTUs_infants, 'Short_sample_ID')

PPVs <- p_vir_frac_infants[["PPB_bacteria"]]

names(PPVs) <- VLP_metadata$Alex_ID_short[match(names(PPVs), VLP_metadata$Individual_ID)]

transmission_p_infant <- transmission_p_infant[ names(transmission_p_infant) %in% names(PPVs) ]
PPVs <- PPVs[ names(PPVs) %in% names(transmission_p_infant) ]

transmission_p_infant_PPV <- transmission_p_infant

for (infant in names(transmission_p_infant)) {
  
  transmission_p_infant_PPV[[infant]] <- transmission_p_infant_PPV[[ infant ]]
  
  transmission_p_infant_PPV[[infant]] <- transmission_p_infant_PPV[[ infant ]][ transmission_p_infant_PPV[[ infant ]] %in% PPVs[[ infant ]] ]
  
}

transmission_p_infant_PPV <- transmission_p_infant_PPV[lapply(transmission_p_infant_PPV,length)>0]

unique(gsub('_.*','',selected_viruses[selected_viruses$Virus %in% unique(unname(unlist(transmission_p_infant_PPV))),]$Host_species))

