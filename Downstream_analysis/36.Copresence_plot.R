setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Copresence plots for lysogenic transmitted viruses
#############################################################

##############################
# Functions
##############################
make_triangles <- function(x, y, point = "up") {
  x <- as.integer(as.factor((x)))
  y <- as.integer(as.factor((y)))
  
  if (point == "up") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x - 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y + 0.5, y + 0.5)
    }, simplify = FALSE)
  } else if (point == "down") {
    newx <- sapply(x, function(x) {
      c(x - 0.5, x + 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y - 0.5, y + 0.5)
    }, simplify = FALSE)
  }
  data.frame(x = unlist(newx), y = unlist(newy))
}

make_squares <- function(x, y) {
  x <- as.integer(as.factor((x)))
  y <- as.integer(as.factor((y)))
  

    newx <- sapply(x, function(x) {
      c(x - 0.5, x - 0.5, x + 0.5, x + 0.5)
    }, simplify = FALSE)
    newy <- sapply(y, function(y) {
      c(y - 0.5, y + 0.5, y + 0.5, y - 0.5)
    }, simplify = FALSE)
  data.frame(x = unlist(newx), y = unlist(newy))
}
##############################
# Loading libraries
##############################
library(tidyverse)
library(patchwork)
##############################
# Input data
##############################
#combo_metadata <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/metadata_combined_for_exp.txt', sep='\t', header=T)

VLP_metadata <- read.table("02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt", sep='\t', header=T)
MGS_metadata <- read.table("02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt", sep='\t', header=T)
common_columns <- c("Universal_fecal_ID", "SAMPLE_ID", "Old_ID", "NG_ID", "NEXT_ID", "Type", "Timepoint", "FAM_ID", "Short_sample_ID", "Short_sample_ID_bact", "Individual_ID")

metadata <- rbind(VLP_metadata[,common_columns], MGS_metadata[,common_columns])
metadata$source <- lapply(metadata$Short_sample_ID, function(x){
                          ifelse(length(grep('V', x))!=0, "VLP", "MGS")
                                                               }   )
metadata$Alex_ID <- paste0(metadata$FAM_ID, '_', metadata$Type, '_', substr(metadata$Short_sample_ID, 1,1), '_', metadata$source, '_', metadata$Timepoint)

selected_viruses <- read.table('02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', header=T)
selected_viruses <- selected_viruses[!is.na(selected_viruses$Virus),]
selected_viruses <- selected_viruses[selected_viruses$temperate==1 & selected_viruses$Transmitted_in_N_related!=0 & !is.na(selected_viruses$Transmitted_in_N_related),]

RPKM_counts_VLP <- read.table('02.CLEAN_DATA/RPKM_counts_VLP.txt', sep='\t', header=T)
RPKM_counts_MGS <- read.table('02.CLEAN_DATA/RPKM_counts_MGS.txt', sep='\t', header=T)

RPKM_combined <- merge(RPKM_counts_VLP, RPKM_counts_MGS, by='row.names', all = T)
row.names(RPKM_combined) <- RPKM_combined$Row.names
RPKM_combined$Row.names <- NULL
RPKM_combined[is.na(RPKM_combined)] <- 0
RPKM_combined <- as.data.frame(t(as.data.frame(t(RPKM_combined))/colSums(RPKM_combined)))*100
RPKM_combined <- RPKM_combined[row.names(RPKM_combined) %in% unique(selected_viruses$Virus),]
row.names(RPKM_combined) <- selected_viruses$ContigID_easy[match(row.names(RPKM_combined), selected_viruses$Virus)]
colnames(RPKM_combined) <- metadata$Alex_ID[match(colnames(RPKM_combined), metadata$Short_sample_ID)]

RPKM_combined_0.95 <- read.table('03.SCRIPTS/NEXT_pilot_FUP_bf_origin/Data_for_Alex/RPKM_counts_combined_0.95_UPD_final_for_exp.txt', sep='\t', header=T)
RPKM_combined_0.95 <- RPKM_combined_0.95[row.names(RPKM_combined_0.95) %in% unique(selected_viruses$Virus),]

row.names(RPKM_combined_0.95) <- selected_viruses$ContigID_easy[match(row.names(RPKM_combined_0.95), selected_viruses$Virus)]
colnames(RPKM_combined_0.95) <- metadata$Alex_ID[match(colnames(RPKM_combined_0.95), metadata$Short_sample_ID)]

RPKM_combined_0.95$MOCK <- 0

microbiome <- read.table("02.CLEAN_DATA/Microbiome_species_unfiltred.txt", sep='\t', header=T)
colnames(microbiome) <- metadata$Alex_ID[match(colnames(microbiome), metadata$Short_sample_ID_bact)]
row.names(microbiome) <- gsub('_', ' ', gsub('.*s__', '', row.names(microbiome)))
microbiome <- microbiome[grep('Bifidobacterium', row.names(microbiome)),]
##############################
# ANALYSIS
##############################

FAM_remove <- list()

for (virusName in row.names(RPKM_combined_0.95)) {
  
  FAM_remove_p_virus <- character()
  
  for (FAM in unique(substr(colnames(RPKM_combined_0.95), 1, 7)) ) {
    
    zerosum_infants <- rowSums( RPKM_combined_0.95[virusName, c(grep(paste0(FAM, '_Infant'), colnames(RPKM_combined_0.95)), ncol(RPKM_combined_0.95), ncol(RPKM_combined_0.95)) ] ) == 0
    zerosum_mothers <- rowSums( RPKM_combined_0.95[virusName, c(grep(paste0(FAM, '_Mother'), colnames(RPKM_combined_0.95)), ncol(RPKM_combined_0.95), ncol(RPKM_combined_0.95)) ] ) == 0
    
    if ( zerosum_infants | zerosum_mothers ) {
      
      FAM_remove_p_virus <- c(FAM_remove_p_virus, FAM)
      
    }
    
  }
  
  FAM_remove[[virusName]] <- FAM_remove_p_virus
  
}


maternal_timepoints <- unique(metadata[metadata$Type=='Mother',]$Timepoint)
infant_timepoints <- c(unique(metadata[metadata$Type=='Infant',]$Timepoint), "B")
sources <- c("VLP", "MGS")
new_column_names <- c("FAM", "Type", "F_member", "Sample_type", "Timepoint")

forplots <- list()

for ( virusName in row.names(RPKM_combined) ) {
  
  virusN <- as.data.frame(t(RPKM_combined[virusName,]))
  
  virusN <- virusN[grep( paste0( FAM_remove[[virusName]], collapse = "|" ), row.names(virusN), invert = T  ), , drop=F]
  
  for ( FAM in unique(substr(row.names(virusN), 1, 7)) ) {
    
    for (sequences in sources ) {
      
      # mothers:
      for (mother_tp in maternal_timepoints) {
        
        if ( length(grep(paste0(FAM, '_Mother_A_', sequences, '_',mother_tp), row.names(virusN)))==0    ) {
          
          virusN[paste0(FAM, '_Mother_A_', sequences, '_',mother_tp),1] <- NA
          
        }
        
      }
      
      # infants:
      for (infant_tp in infant_timepoints) {
        
        if ( length(grep(paste0(FAM, '_Infant_C_', sequences, '_', infant_tp, '$'), row.names(virusN)))==0    ) {
          
          virusN[paste0(FAM, '_Infant_C_', sequences, '_', infant_tp), 1] <- NA
          
        }
        # checking if there is a twin
        
        if ( length( grep( paste0(FAM, '_Infant_D_'),  row.names(virusN)) ) != 0  ) {
          
          if ( length(grep(paste0(FAM, '_Infant_D_', sequences, '_', infant_tp, '$'), row.names(virusN)))==0    ) {
            
            virusN[paste0(FAM, '_Infant_D_', sequences, '_', infant_tp), 1] <- NA
            
          }
          
        }
        
      }
      
    }
    
  }
  
  virusN$names <- row.names(virusN)
  
  virusN <- separate(virusN, names, into = new_column_names, sep = "_")
  
  colnames(virusN)[1] <- "Presence"
  
  virusN[!is.na(virusN$Presence) & virusN$Presence > 0 ,]$Presence <- 1
  
  virusN$Type <- paste0(virusN$Type, '_', virusN$F_member)
  virusN$F_member <- NULL
  
  virusN_wide <- virusN %>% pivot_wider(names_from = "Sample_type", values_from = "Presence")
  
  virusN_wide$Timepoint <- droplevels(factor(virusN_wide$Timepoint, levels = c('P3', 'P7', 'B', 'M1', 'M2', 'M3', 'M6', 'M9', 'M12')))
  
  virusN_wide$Type <- paste0(virusN_wide$FAM, '_', virusN_wide$Type)
  
  newcoord_up <- make_triangles(virusN_wide$Timepoint, virusN_wide$Type)
  newcoord_down <- make_triangles(virusN_wide$Timepoint, virusN_wide$Type, point = "down")
  newcoord_down <- newcoord_down %>% select(xdown = x, ydown = y)
  
  repdata <- map_df(1:nrow(virusN_wide), function(i) virusN_wide[rep(i, 3), ])
  newdata <- bind_cols(repdata, newcoord_up, newcoord_down)
  newdata$Type <- factor(newdata$Type, ordered = T)
  newdata$MGS <- as.character(newdata$MGS)
  newdata$VLP <- as.character(newdata$VLP)
  
  for (FAM in sort(unique(newdata$FAM)) ) {
    
    newdata[newdata$FAM==FAM, "y"] <- newdata[newdata$FAM==FAM, "y"] + as.numeric(grep(FAM, sort(unique(newdata$FAM)))) - 0.5
    newdata[newdata$FAM==FAM, "ydown"] <- newdata[newdata$FAM==FAM, "ydown"] + as.numeric(grep(FAM, sort(unique(newdata$FAM)))) - 0.5
  
    }
  
  newdata$Breaks <- NA
  
  for (i in unique(newdata$Type)) {
    
    A <- as.numeric(unname(unlist(unique(newdata[newdata$Type==i,"y"]))))
    
    newdata[newdata$Type==i,"Breaks"] <- (A[-1] + A[-length(A)]) / 2
    
  }
  
  #newdata[,c(6:10)] <- newdata[,c(6:10)]/2
  forplots[[virusName]] <- newdata
  
}

forbacplots <- list()

for ( bacteriumName in row.names(microbiome) ) {
  
  bacteriumN <- as.data.frame(t(microbiome[bacteriumName,]))
  
  for ( FAM in unique(substr(row.names(bacteriumN), 1, 7)) ) {
    
    #for (sequences in sources ) {
      
      # mothers:
      for (mother_tp in maternal_timepoints) {
        
        if ( length(grep(paste0(FAM, '_Mother_A_MGS_',mother_tp), row.names(bacteriumN)))==0    ) {
          
          bacteriumN[paste0(FAM, '_Mother_A_MGS_',mother_tp),1] <- NA
          
        }
        
      }
      
      # infants:
      for (infant_tp in infant_timepoints) {
        
        if ( length(grep(paste0(FAM, '_Infant_C_MGS_', infant_tp, '$'), row.names(bacteriumN)))==0    ) {
          
          bacteriumN[paste0(FAM, '_Infant_C_MGS_', infant_tp), 1] <- NA
          
        }
        # checking if there is a twin
        
        if ( length( grep( paste0(FAM, '_Infant_D_'),  row.names(bacteriumN)) ) != 0  ) {
          
          if ( length(grep(paste0(FAM, '_Infant_D_MGS_', infant_tp, '$'), row.names(bacteriumN)))==0    ) {
            
            bacteriumN[paste0(FAM, '_Infant_D_MGS_', infant_tp), 1] <- NA
            
          }
          
        }
        
      }
      
    #}
    
  }
  
  bacteriumN$names <- row.names(bacteriumN)
  
  bacteriumN <- separate(bacteriumN, names, into = new_column_names, sep = "_")
  
  colnames(bacteriumN)[1] <- "Presence"
  
  bacteriumN[!is.na(bacteriumN$Presence) & bacteriumN$Presence > 0 ,]$Presence <- 1
  
  bacteriumN$Type <- paste0(bacteriumN$Type, '_', bacteriumN$F_member)
  bacteriumN$F_member <- NULL
  
  bacteriumN_wide <- bacteriumN %>% pivot_wider(names_from = "Sample_type", values_from = "Presence")
  
  bacteriumN_wide$Timepoint <- droplevels(factor(bacteriumN_wide$Timepoint, levels = c('P3', 'P7', 'B', 'M1', 'M2', 'M3', 'M6', 'M9', 'M12')))
  
  bacteriumN_wide$Type <- paste0(bacteriumN_wide$FAM, '_', bacteriumN_wide$Type)
  
  newcoord <- make_squares(bacteriumN_wide$Timepoint, bacteriumN_wide$Type)

  repdata <- map_df(1:nrow(bacteriumN_wide), function(i) bacteriumN_wide[rep(i, 4), ])
  newdata <- bind_cols(repdata, newcoord)
  newdata$Type <- factor(newdata$Type, ordered = T)
  newdata$MGS <- as.character(newdata$MGS)

  for (FAM in sort(unique(newdata$FAM)) ) {
    
    newdata[newdata$FAM==FAM, "y"] <- newdata[newdata$FAM==FAM, "y"] + as.numeric(grep(FAM, sort(unique(newdata$FAM)))) - 0.5
    
  }
  
  newdata$Breaks <- NA
  
  for (i in unique(newdata$Type)) {
    
    A <- as.numeric(unname(unlist(unique(newdata[newdata$Type==i,"y"]))))
    
    newdata[newdata$Type==i,"Breaks"] <- (A[-1] + A[-length(A)]) / 2
    
  }
  
  #newdata[,c(6:10)] <- newdata[,c(6:10)]/2
  forbacplots[[bacteriumName]] <- newdata
  
}

tmp0 <- forbacplots[[4]][grep('FAM0234', forbacplots[[4]]$Type),]
tmp0$VLP <- NA
tmp0$xdown <- NA
tmp0$ydown <- NA
tmp0 <- tmp0[,c(1:3,8,4:6,9,10,7)]
tmp0[tmp0$Type=='FAM0234_Mother_A',]$y <- tmp0[tmp0$Type=='FAM0234_Mother_A',]$y - 38
tmp0[tmp0$Type=='FAM0234_Mother_A',]$Breaks <- tmp0[tmp0$Type=='FAM0234_Mother_A',]$Breaks - 38
tmp0[tmp0$Type=='FAM0234_Infant_C',]$y <- tmp0[tmp0$Type=='FAM0234_Infant_C',]$y - 40
tmp0[tmp0$Type=='FAM0234_Infant_C',]$Breaks <- tmp0[tmp0$Type=='FAM0234_Infant_C',]$Breaks - 40
tmp0$Type <- as.character(tmp0$Type)
tmp0[tmp0$Type=='FAM0234_Mother_A',]$Type <- 'FAM0234_Mother_B'
tmp0[tmp0$Type=='FAM0234_Infant_C',]$Type <- 'FAM0234_Infant_B'


pictograms_all <- list()

for ( virusName in names(forplots) ) {
  
  pictograms_all[[virusName]] <- ggplot(forplots[[virusName]]) +
    geom_polygon(aes(x = x, y = y, fill = MGS, group = interaction(Type, Timepoint)), color = "black") +
    scale_fill_manual(values=c('white', '#D72323', '#7F7F7F'),  limits = c("0", "1", "NA")) +
    #ggnewscale::new_scale_fill() +
    labs(x="Timepoint", y="Families") +
    geom_polygon(aes(x = xdown, y = ydown, fill = VLP, group = interaction(Type, Timepoint)), color = "black") +
    scale_fill_manual(values=c('white', '#D72323', '#7F7F7F'),  limits = c("0", "1", "NA")) +
    ggtitle(names(forplots[virusName])) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                       labels = levels(forplots[[virusName]]$Timepoint)) +
    scale_y_continuous(breaks = sort(unique(forplots[[virusName]]$Breaks)),
                       labels = levels(forplots[[virusName]]$Type)) +
    coord_equal() + 
    theme_linedraw() +
    theme(panel.grid = element_blank(), 
          axis.title = element_text(face="bold")) + 
    guides(fill= guide_legend("Presence", title.position = "top", title.hjust = 0.5))
  
}


combined_plot <- (pictograms_all[[1]] / pictograms_all[[2]] / pictograms_all[[5]]) | (pictograms_all[[3]] / pictograms_all[[4]] )


combined_plot <- combined_plot +
  plot_layout(ncol=2, guides = "collect") + 
  plot_annotation(title = "") & theme(legend.position = "bottom", 
                                      axis.title = element_text(size=6),
                                      axis.text = element_text(size=5), 
                                      plot.title = element_text(size=5),
                                      legend.title = element_text(size=7),
                                      legend.text = element_text(size=7)) 



tmp <- rbind(forplots[["L34922_LS1"]], tmp0)


ggplot(tmp) +
  geom_polygon(aes(x = x, y = y, fill = MGS, group = interaction(Type, Timepoint)), color = "black", show.legend = F) +
  scale_fill_manual(values=c('white', '#D72323', '#7F7F7F'),  limits = c("0", "1", "NA")) +
  ggnewscale::new_scale_fill() +
  labs(x="Timepoint", y="") +
  geom_polygon(aes(x = xdown, y = ydown, fill = VLP, group = interaction(Type, Timepoint)), color = "black") +
  scale_fill_manual(values=c('white', '#D72323', '#7F7F7F'),  limits = c("0", "1", "NA")) +
  ggtitle("Prsenece patterns") +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9), 
                     labels = levels(tmp$Timepoint)) +
  scale_y_continuous(breaks = c(1,3),
                     labels = c("Infant","Mother")) +
  coord_equal() + 
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(face="bold"),
        legend.position = "bottom") + 
  guides(fill= guide_legend("Presence", title.position = "top", title.hjust = 0.5))






# pdf('04.PLOTS/Copresence_transmitted_lysogenoc.pdf', width=15/2.54, height=13/2.54)
# combined_plot
# dev.off()

##### FOR VISUALIZATION ####
write.table(tmp, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig6B_L34922_LS1_copresence.txt", sep='\t', row.names = F)
