setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Here we explore co-transmission of viruses and their 
# bacterial hosts
#############################################################

##############################
# Functions
##############################
flatten_correlation_matrix <- function(cor_matrix, pvalue_matrix) {
  # Get the row names and column names
  row_names <- row.names(cor_matrix)
  col_names <- colnames(cor_matrix)
  
  # Initialize an empty data frame to store the flattened matrix
  flattened <- data.frame(row = character(),
                          col = character(),
                          correlation = numeric(),
                          pvalue = numeric(),
                          stringsAsFactors = FALSE)
  
  # Iterate over each cell in the correlation matrix
  for (i in seq_along(row_names)) {
    for (j in seq_along(col_names)) {
      # Get the row name, column name, correlation value, and p-value
      row_name <- row_names[i]
      col_name <- col_names[j]
      correlation <- cor_matrix[i, j]
      pvalue <- pvalue_matrix[i, j]
      
      # Check if the correlation value and p-value are not NA
      if (!is.na(correlation) && !is.na(pvalue)) {
        # Create a new row in the flattened data frame
        new_row <- data.frame(row = row_name, col = col_name,
                              correlation = correlation, pvalue = pvalue)
        
        # Append the new row to the flattened data frame
        flattened <- rbind(flattened, new_row)
      }
    }
  }
  
  # Return the flattened data frame
  return(flattened)
}
##############################
# Loading libraries
##############################
library(vegan)
library(corrplot)
library(reshape2)
library(ggmosaic)
library(tidyverse)
##############################
# Input data
##############################

selected_viruses <- read.table("02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt", sep='\t', header=T)

# selecting viruses and their hosts for the experiment:
check_cotransmission <- selected_viruses
# filter 0: select only those viruses that were used for the strain sharing analysis: 
check_cotransmission <- check_cotransmission[!is.na(check_cotransmission$Distances_Related_lower) & check_cotransmission$Distances_Related_lower=="YES",]
# filter 1: do not consider those viruses for which at least 3 transmission/strain sharing events were detected
check_cotransmission <- check_cotransmission[check_cotransmission$Transmitted_in_N_related > 2,]
# filter 2: select only those hosts that were used for the strain sharing analysis: 
check_cotransmission <- check_cotransmission[!is.na(check_cotransmission$Host_Distances_Related_lower) & check_cotransmission$Host_Distances_Related_lower=="YES", ]
# filter 3: select only those hosts for which at least 3 transmission/strain sharing events were detected
check_cotransmission <- check_cotransmission[check_cotransmission$Host_Transmitted_in_N_related > 2, ]


# cutpoint: 
check_cotransmission$cutpoint_virus <- ifelse(check_cotransmission$Youden_index >= check_cotransmission$FDR_ipv_Youden, check_cotransmission$FDR_ipv_Youden, check_cotransmission$Youden_index)
check_cotransmission$cutpoint_bacterium <- ifelse(check_cotransmission$Host_Youden_index >= check_cotransmission$Host_FDR_ipv_Youden, check_cotransmission$Host_FDR_ipv_Youden, check_cotransmission$Host_Youden_index)

# alternative easy name for host:
check_cotransmission$Host_easy <- paste0(check_cotransmission$Host_species, '_', gsub('SGB', '', check_cotransmission$Host_SGB) )

check_cotransmission <- check_cotransmission[order(check_cotransmission$Host_easy),]

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)
virus <- virus[names(virus) %in% check_cotransmission$Virus]

#### normalizing distance matrices by the median, as thresholds are calculated using the normalized distances:
virus <- lapply(virus, function(x) {
  x <- x/median( x[upper.tri(x)] )
})

#### exchanging distances values with 0 and 1 depending on the threshold for strain sharing event:
for (virusName in names(virus)) {
  
  virusN <- virus[[virusName]]
  
  virusN[virusN <= check_cotransmission[check_cotransmission$Virus==virusName,]$cutpoint_virus] <- 0
  
  virusN[virusN > check_cotransmission[check_cotransmission$Virus==virusName,]$cutpoint_virus] <- 1
  
  virus[[virusName]] <- virusN
}

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_', MGS_metadata$source, '_', MGS_metadata$Timepoint)

folder = "01.RAW_DATA/strainphlan_4_distance_matrix/"  
file_list = list.files(path=folder, pattern="*.csv")  
bacterium=lapply(paste0(folder, file_list), function(x) read.csv(x, header=T))
names(bacterium) <- gsub("_dmat_.csv", '', file_list)
bacterium <- bacterium[names(bacterium) %in% check_cotransmission$Host_SGB]
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
  
  bacteriumN[bacteriumN <= check_cotransmission[check_cotransmission$Host_SGB==bacteriumName,]$cutpoint_bacterium] <- 0
  
  bacteriumN[bacteriumN > check_cotransmission[check_cotransmission$Host_SGB==bacteriumName,]$cutpoint_bacterium] <- 1
  
  bacterium[[bacteriumName]] <- bacteriumN
}

set.seed(888)
##############################
# ANALYSIS
##############################

co_transmission_cor <- data.frame(matrix(NA, ncol=length(unique(check_cotransmission$Virus)),
                                     nrow=length(unique(check_cotransmission$Host_SGB)))  )

colnames(co_transmission_cor) <- unique(check_cotransmission$Virus)
row.names(co_transmission_cor) <- unique(check_cotransmission$Host_SGB)

co_transmission_pval <- co_transmission_cor

mantel_partial_results <- list()

# calculating the correlation between matrices of strain sharing reconstructed for viruses and their predicted hosts:

# Iterate over each virus in the dataset
for (virusName in unique(check_cotransmission$Virus)) {
  
  # Obtain the virus genetic distance matrix
  virusN <- virus[[virusName]]
  
  # Parse sample names to handle cases where virus strain was reconstructed in both MGS and VLP
  sample_names <- data.frame(colnames(virusN))
  colnames(sample_names) <- 'sample'
  sample_names$Var1 <- gsub("VLP_|MGS_", "", sample_names$sample)
  sample_names$source <- NA
  sample_names[grep('VLP', sample_names$sample),]$source <- 'VLP'
  sample_names[grep('MGS', sample_names$sample),]$source <- 'MGS'
  repeated <- data.frame(table(gsub("VLP_|MGS_", "", colnames(virusN))))
  sample_names <- merge(sample_names, repeated, by='Var1')
  sample_names <- sample_names[sample_names$Freq==1 | sample_names$source=='VLP',]
  sample_names$bacterial <- gsub('VLP', 'MGS', sample_names$sample)
  
  # Iterate over each bacterial host strain predicted to be infected by the virus:
  for (bacteriumName in unique(check_cotransmission$Host_SGB)) {
    
    # Obtain the bacterial distance matrix
    bacteriumN <- bacterium[[bacteriumName]]
    
    # Merge names of samples with reconstructed virus and bacterial strains
    sample_names_bacterial <- data.frame(colnames(bacteriumN))
    colnames(sample_names_bacterial) <- 'bacterial'
    sample_names_concurrent <- merge(sample_names, sample_names_bacterial, by='bacterial')
    
    # Extract overlapping samples where both virus and bacterial strains were reconstructed
    virusN_concurrent <- virusN[sample_names_concurrent$sample, sample_names_concurrent$sample]
    bacteriumN <- bacteriumN[sample_names_concurrent$bacterial, sample_names_concurrent$bacterial]
    
    # the number of concurrent samples should be larger than 7, even though 7! is larger than 1000
    if (length(sample_names_concurrent$bacterial)>7) {
      
      dist_matrix_contol <- virusN_concurrent
      
      for (i in colnames(dist_matrix_contol)) {
        
        for (j in row.names(dist_matrix_contol)) {
          
          if ( substr(i, 1, 16)==substr(j, 1, 16) ) {
            dist_matrix_contol[i,j] <- 0
          } else {
            dist_matrix_contol[i,j] <- 1
          }
          
        }
        
      }
      
      dist_matrix_contol_NA <- dist_matrix_contol
      dist_matrix_contol_NA[is.na(bacteriumN)] <- NA
      
      if ( setequal(virusN_concurrent, bacteriumN) == T ) { # if topology is the same
        
        mantel_partial_results[[paste0(virusName, '_', bacteriumName)]] <- mantel(virusN_concurrent, dist_matrix_contol, method = "pearson", permutations = 999, na.rm = T)
        
        co_transmission_cor[bacteriumName,virusName] <- 1
        co_transmission_pval[bacteriumName,virusName] <- mantel_partial_results[[paste0(virusName, '_', bacteriumName)]]$signif
         
      } else {
        
        if (setequal(bacteriumN, dist_matrix_contol) == T | setequal(virusN_concurrent, dist_matrix_contol) | setequal(dist_matrix_contol_NA, bacteriumN) ) {
          
          mantel_partial_results[[paste0(virusName, '_', bacteriumName)]] <- "No time adjustment possible"
          
        } else {
          
          mantel_partial_results[[paste0(virusName, '_', bacteriumName)]] <- mantel.partial(virusN_concurrent, bacteriumN, dist_matrix_contol, method = "pearson", permutations = 999, na.rm = T)
          
          # Store the correlation coefficient and p-value
          co_transmission_cor[bacteriumName,virusName] <- mantel_partial_results[[paste0(virusName, '_', bacteriumName)]]$statistic
          co_transmission_pval[bacteriumName,virusName] <- mantel_partial_results[[paste0(virusName, '_', bacteriumName)]]$signif
          
        }
        
        
      }

      
    } else {
      mantel_partial_results[[paste0(virusName, '_', bacteriumName)]] <- 'No concurrent samples'
    }

  }
  
}

colnames(co_transmission_cor) <- check_cotransmission$ContigID_easy[match(colnames(co_transmission_cor), check_cotransmission$Virus)]
row.names(co_transmission_cor) <- check_cotransmission$Host_easy[match(row.names(co_transmission_cor), check_cotransmission$Host_SGB)]

colnames(co_transmission_pval) <- check_cotransmission$ContigID_easy[match(colnames(co_transmission_pval), check_cotransmission$Virus)]
row.names(co_transmission_pval) <- check_cotransmission$Host_easy[match(row.names(co_transmission_pval), check_cotransmission$Host_SGB)]

co_transmission_cor <- co_transmission_cor[,colSums(is.na(co_transmission_cor)) != nrow(co_transmission_cor)]
co_transmission_cor <- co_transmission_cor[rowSums(is.na(co_transmission_cor))!=ncol(co_transmission_cor),]

co_transmission_cor <- co_transmission_cor[order(row.names(co_transmission_cor)),]
co_transmission_pval <- co_transmission_pval[,colnames(co_transmission_pval) %in% colnames(co_transmission_cor)]
co_transmission_pval <- co_transmission_pval[order(row.names(co_transmission_pval)),]

co_transmission_cor <- as.matrix(co_transmission_cor)
co_transmission_pval <- as.matrix(co_transmission_pval)

# multiple test correction:
cotransmission_flat <- flatten_correlation_matrix(co_transmission_cor, co_transmission_pval)
cotransmission_flat$FDR <- p.adjust(cotransmission_flat$pvalue, method = "BH")
co_transmission_FDR <- co_transmission_pval

for (i in 1:nrow(cotransmission_flat)) {
  
  virusName <- cotransmission_flat[ i,  'col']
  bacteriumName <- cotransmission_flat[ i,  'row']
  
  co_transmission_FDR[bacteriumName,virusName] <- cotransmission_flat[i, 'FDR']
  
}

co_transmission_cor[is.nan(co_transmission_cor) | !is.finite(co_transmission_cor)] <- NA 
co_transmission_FDR[is.nan(co_transmission_FDR) | !is.finite(co_transmission_FDR)] <- NA
co_transmission_FDR[is.na(co_transmission_FDR)] <- NA

bgcolors <- matrix(0, nrow = nrow(co_transmission_cor), 
                   ncol = ncol(co_transmission_cor), 
                   dimnames = list(row.names(co_transmission_cor), colnames(co_transmission_cor)))


for (virusName in colnames(co_transmission_cor)) {
  
  bacteriumName <- check_cotransmission[check_cotransmission$ContigID_easy==virusName,]$Host_easy
  
  bgcolors[bacteriumName, virusName] <- NA
}

# some mantel exact stat artifact I guess:
co_transmission_cor[co_transmission_cor>1] <- 1


pdf('./04.PLOTS/Virus_bacteria_co_transmission_corrplot_FDR_0.1.pdf', width=25/2.54, height=22/2.54)
co_transmission_cor[co_transmission_FDR > 0.10] <- 0
co_transmission_cor[is.na(co_transmission_FDR)] <- NA
corrplot(bgcolors, na.label = "square", na.label.col = "#A9C52F", tl.col = "white", cl.pos = "none")
corrplot(co_transmission_cor,
         bg=NA,
         p.mat=co_transmission_FDR, 
         cl.lim=c(0,1), is.corr = F,
         insig = "blank",
         na.label = "X", 
         na.label.col = "lightgrey", 
         tl.col='black', 
         add = T,
         col=colorRampPalette(c("white","lightblue","navy"))(100), cl.pos="b")

dev.off()


# test if co-transmission is enriched in predicted virus-host pairs vs random virus-bacteria pairs:

cotransmission_flat$virus_host_pair <- "NO"
for (i in 1:nrow(cotransmission_flat)) {
  
  bacteriumName <- cotransmission_flat[i,"row"]
  virusName <- cotransmission_flat[i,"col"]
  
  if (bacteriumName %in% check_cotransmission[check_cotransmission$ContigID_easy==virusName,"Host_easy"] ) {
    cotransmission_flat[i, "virus_host_pair"] <- "YES"
  }
  
}

dim(cotransmission_flat[cotransmission_flat$FDR <= 0.05 & cotransmission_flat$virus_host_pair=="YES",])[1]
dim(cotransmission_flat[cotransmission_flat$FDR > 0.05 & cotransmission_flat$virus_host_pair=="YES",])[1]
dim(cotransmission_flat[cotransmission_flat$FDR <= 0.05 & cotransmission_flat$virus_host_pair=="NO",])[1]
dim(cotransmission_flat[cotransmission_flat$FDR > 0.05 & cotransmission_flat$virus_host_pair=="NO",])[1]

cotransmission_stat <- matrix( c(20, # pair & cotransmitted
                                 25, # pair & not cotransmitted
                                 85, # not a pair & cotransmitted
                                 242), # not a pair & not cotransmittd
                                 nrow=2,
                                 dimnames = list(c("Cotransmitted", "Not-cotransmitted"),
                                             c("Paired", "Unpaired")))


fisher.test(cotransmission_stat, alternative="greater")

#### FOR SUPPLEMENTARY #####
write.table(cotransmission_flat, '05.MANUSCRIPT/Supplementary_tables/Virus_bacteria_cotransmission_stat.txt', sep='\t', quote=F)



cotransmission_stat_melt <- melt(cotransmission_stat)
cotransmission_stat_melt$prop <- c(20/45, 25/45, 85/327, 242/327)
cotransmission_stat_melt$paired.count <- c(45, 45, 327, 327)

ggplot(cotransmission_stat_melt, aes(x=Var2, y=prop, width=paired.count, fill=Var1)) + 
  labs(x="Virus - predicted host pair", y="Proportion", fill="Transmission pattern") + 
  geom_bar(stat = "identity", position = "fill", colour = "black") + 
  geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5)) +
  facet_grid(~Var2, scales = "free_x", space = "free_x")  +
  theme(panel.spacing.x = unit(0, "npc")) +
  scale_fill_brewer(palette = "RdYlGn") +
  theme_linedraw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_text(face="bold"),
        strip.background = element_blank(),
        strip.text.x = element_blank())


##############################
# FOR VISUALIZATION
##############################
write.table(bgcolors, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5A_bg_colors_bacteria_virus_pair.txt", sep='\t', quote=F)
write.table(co_transmission_cor, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5A_cotransmission_correlation_matrix.txt", sep='\t', quote=F)
write.table(co_transmission_FDR, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5A_cotransmission_FDR_matrix.txt", sep='\t', quote=F)
write.table(cotransmission_stat_melt, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/5B_cotransmission_stat.txt", sep='\t', quote=F, row.names=F)


#### Linkage of Disequilibrium

co_transmission_LD <- data.frame(matrix(NA, ncol=length(unique(check_cotransmission$Virus)),
                                         nrow=length(unique(check_cotransmission$Host_SGB)))  )

colnames(co_transmission_LD) <- unique(check_cotransmission$Virus)
row.names(co_transmission_LD) <- unique(check_cotransmission$Host_SGB)

co_transmission_LDpval <- co_transmission_LD
for (virusName in unique(check_cotransmission$Virus)) {
  
  # Obtain the virus genetic distance matrix
  virusN <- virus[[virusName]]
  
  # Parse sample names to handle cases where virus strain was reconstructed in both MGS and VLP
  sample_names <- data.frame(colnames(virusN))
  colnames(sample_names) <- 'sample'
  sample_names$Var1 <- gsub("VLP_|MGS_", "", sample_names$sample)
  sample_names$source <- NA
  sample_names[grep('VLP', sample_names$sample),]$source <- 'VLP'
  sample_names[grep('MGS', sample_names$sample),]$source <- 'MGS'
  repeated <- data.frame(table(gsub("VLP_|MGS_", "", colnames(virusN))))
  sample_names <- merge(sample_names, repeated, by='Var1')
  sample_names <- sample_names[sample_names$Freq==1 | sample_names$source=='VLP',]
  sample_names$bacterial <- gsub('VLP', 'MGS', sample_names$sample)
  
  # Iterate over each bacterial host strain predicted to be infected by the virus:
  for (bacteriumName in unique(check_cotransmission$Host_SGB)) {
    
    # Obtain the bacterial distance matrix
    bacteriumN <- bacterium[[bacteriumName]]
    
    # Merge names of samples with reconstructed virus and bacterial strains
    sample_names_bacterial <- data.frame(colnames(bacteriumN))
    colnames(sample_names_bacterial) <- 'bacterial'
    sample_names_concurrent <- merge(sample_names, sample_names_bacterial, by='bacterial')
    
    # Extract overlapping samples where both virus and bacterial strains were reconstructed
    virusN_concurrent <- virusN[sample_names_concurrent$sample, sample_names_concurrent$sample, drop=F]
    bacteriumN <- bacteriumN[sample_names_concurrent$bacterial, sample_names_concurrent$bacterial, drop=F]
    
    if (ncol(bacteriumN) > 7) {
      A <- sum(virusN_concurrent[upper.tri(virusN_concurrent)]==0) / length(virusN_concurrent[upper.tri(virusN_concurrent)])
      B <- sum(bacteriumN[upper.tri(bacteriumN)]==0) / length((bacteriumN[upper.tri(bacteriumN)]))
      
      num_zeros <- sum(virusN_concurrent[upper.tri(virusN_concurrent)] == 0 & bacteriumN[upper.tri(bacteriumN)] == 0, na.rm = TRUE)
      
      coshared_freq_real <- num_zeros/length(virusN_concurrent[upper.tri(virusN_concurrent)])
      
      coshared_freq_est <- A*B
      
      co_transmission_LD[bacteriumName,virusName] <- coshared_freq_real - coshared_freq_est
      
      chisq_D <- (ncol(virusN_concurrent)*(co_transmission_LD[bacteriumName,virusName])^2) / ((A * (1-A)) * B * (1 - B))
      
      co_transmission_LDpval[bacteriumName,virusName] <- 1 - pchisq(chisq_D, df=1)
      
    }
    
  }
  
}

colnames(co_transmission_LD) <- check_cotransmission$ContigID_easy[match(colnames(co_transmission_LD), check_cotransmission$Virus)]
row.names(co_transmission_LD) <- check_cotransmission$Host_easy[match(row.names(co_transmission_LD), check_cotransmission$Host_SGB)]

colnames(co_transmission_LDpval) <- check_cotransmission$ContigID_easy[match(colnames(co_transmission_LDpval), check_cotransmission$Virus)]
row.names(co_transmission_LDpval) <- check_cotransmission$Host_easy[match(row.names(co_transmission_LDpval), check_cotransmission$Host_SGB)]

#co_transmission_LD <- co_transmission_LD[,colSums(is.na(co_transmission_LD)) != nrow(co_transmission_LD)]
#co_transmission_LD <- co_transmission_LD[rowSums(is.na(co_transmission_LD)) != ncol(co_transmission_LD),]

#co_transmission_LD <- co_transmission_LD[order(row.names(co_transmission_LD)),]
#co_transmission_LDpval <- co_transmission_LDpval[,colnames(co_transmission_LDpval) %in% colnames(co_transmission_LD)]
#co_transmission_LDpval <- co_transmission_LDpval[order(row.names(co_transmission_LDpval)),]

co_transmission_LD <- as.matrix(co_transmission_LD)
co_transmission_LDpval <- as.matrix(co_transmission_LDpval)

# multiple test correction:
cotransmissionLD_flat <- flatten_correlation_matrix(co_transmission_LD, co_transmission_LDpval)
cotransmissionLD_flat$FDR <- p.adjust(cotransmissionLD_flat$pvalue, method = "BH")
co_transmissionLD_FDR <- co_transmission_LDpval

for (i in 1:nrow(cotransmissionLD_flat)) {
  
  virusName <- cotransmissionLD_flat[ i,  'col']
  bacteriumName <- cotransmissionLD_flat[ i,  'row']
  
  co_transmissionLD_FDR[bacteriumName,virusName] <- cotransmissionLD_flat[i, 'FDR']
  
}

colnames(cotransmissionLD_flat)[1] <- 'Bacterium'
colnames(cotransmissionLD_flat)[2] <- 'Virus'
colnames(cotransmissionLD_flat)[3] <- 'LD'

#### FOR SUPPLEMENTARY #####
write.table(cotransmissionLD_flat, '05.MANUSCRIPT/Supplementary_tables/Virus_bacteria_cotransmission_LD_stat.txt', sep='\t', quote=F, row.names = F)


for (virusName in colnames(co_transmission_LD)) {
  
  bacteriumName <- check_cotransmission[check_cotransmission$ContigID_easy==virusName,]$Host_easy
  
  bgcolors[bacteriumName, virusName] <- NA
}


co_transmission_LD[co_transmissionLD_FDR > 0.05] <- 0
co_transmission_LD[is.na(co_transmissionLD_FDR)] <- NA
corrplot(bgcolors, na.label = "square", na.label.col = "#A9C52F", tl.col = "white", cl.pos = "none")
corrplot(co_transmission_LD,
         bg=NA,
         p.mat=co_transmission_FDR, 
         cl.lim=c(0,1), is.corr = F,
         insig = "blank",
         na.label = "X", 
         na.label.col = "lightgrey", 
         tl.col='black', 
         add = T,
         col=colorRampPalette(c("white","lightblue","navy"))(100), cl.pos="b")


cotransmissionLD_flat$virus_host_pair <- "NO"
for (i in 1:nrow(cotransmissionLD_flat)) {
  
  bacteriumName <- cotransmissionLD_flat[i,"row"]
  virusName <- cotransmissionLD_flat[i,"col"]
  
  if (bacteriumName %in% check_cotransmission[check_cotransmission$ContigID_easy==virusName,"Host_easy"] ) {
    cotransmissionLD_flat[i, "virus_host_pair"] <- "YES"
  }
  
}

cotransmissionLD_flat$Significane_level <- NA
cotransmissionLD_flat[cotransmissionLD_flat$FDR > 0.05,]$Significane_level <- ""
cotransmissionLD_flat[cotransmissionLD_flat$FDR <= 0.05,]$Significane_level <- "*"
cotransmissionLD_flat[cotransmissionLD_flat$FDR <= 0.01,]$Significane_level <- "**"
cotransmissionLD_flat[cotransmissionLD_flat$FDR <= 0.001,]$Significane_level <- "***"


ggplot(cotransmissionLD_flat, aes(row, col, fill = correlation))+
  labs(x='Bacterial strains', y='Virus strain') +
  
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.25,0.25), space = "Lab", 
                       name="LD") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text( vjust = 1, 
                                    size = 10, hjust = 1)) +
  
  geom_tile(aes(color=virus_host_pair), width = 0.98, height = 0.98) +
  scale_color_manual(guide = FALSE, values = c("YES" = "black", "NO" = "white")) +
  geom_text(aes(label=Significane_level), color="black", size=3) + 
  coord_flip()

boxplot(cotransmissionLD_flat$correlation ~ cotransmissionLD_flat$virus_host_pair)

dim(cotransmissionLD_flat[cotransmissionLD_flat$FDR <= 0.05 & cotransmissionLD_flat$virus_host_pair=="YES",])[1]
dim(cotransmissionLD_flat[cotransmissionLD_flat$FDR > 0.05 & cotransmissionLD_flat$virus_host_pair=="YES",])[1]
dim(cotransmissionLD_flat[cotransmissionLD_flat$FDR <= 0.05 & cotransmissionLD_flat$virus_host_pair=="NO",])[1]
dim(cotransmissionLD_flat[cotransmissionLD_flat$FDR > 0.05 & cotransmissionLD_flat$virus_host_pair=="NO",])[1]

cotransmissionLD_stat <- matrix( c(14, # pair & cotransmitted
                                 18, # pair & not cotransmitted
                                 65, # not a pair & cotransmitted
                                 198), # not a pair & not cotransmittd
                               nrow=2,
                               dimnames = list(c("Cotransmitted", "Not-cotransmitted"),
                                               c("Paired", "Unpaired")))

fisher.test(cotransmissionLD_stat, alternative="greater")

mosaicplot(cotransmissionLD_stat, color=T)
#### calculating the background: are bacterial strains co-transmitted with each other? 

 cotransmission_BB_cor <- matrix(NA, ncol=length(unique(check_cotransmission$Host_SGB)),
                                              nrow=length(unique(check_cotransmission$Host_SGB)),
                                 dimnames = list(unique(check_cotransmission$Host_SGB), unique(check_cotransmission$Host_SGB)))  

 cotransmission_BB_pval <-  cotransmission_BB_cor
 
 mantel_partial_BB_results <- list()
 
 for (bacteriumName in colnames(cotransmission_BB_cor)) {
   
   # Obtain the virus genetic distance matrix
   bacteriumN <- bacterium[[bacteriumName]]
   
   # Parse sample names to handle cases where virus strain was reconstructed in both MGS and VLP
   sample_names <- colnames(bacteriumN)
   
   # Iterate over each bacterial host strain predicted to be infected by the virus:
   for (bacteriumName2 in colnames(cotransmission_BB_cor)) {
     
     # Obtain the bacterial distance matrix
     bacteriumN2 <- bacterium[[bacteriumName2]]
     
     # Merge names of samples with reconstructed virus and bacterial strains
     sample_names_bacterial <- colnames(bacteriumN2)
     
     sample_names_concurrent <- intersect(sample_names, sample_names_bacterial)
     
     # Extract overlapping samples where both virus and bacterial strains were reconstructed
     bacteriumN_concurrent <- bacteriumN[sample_names_concurrent, sample_names_concurrent]
     bacteriumN2 <- bacteriumN2[sample_names_concurrent, sample_names_concurrent]
     
     if (length(sample_names_concurrent)>=7) {
       
       dist_matrix_contol <- bacteriumN_concurrent
       
       for (i in colnames(dist_matrix_contol)) {
         
         for (j in row.names(dist_matrix_contol)) {
           
           if ( substr(i, 1, 16)==substr(j, 1, 16) ) {
             dist_matrix_contol[i,j] <- 0
           } else {
             dist_matrix_contol[i,j] <- 1
           }
           
         }
         
       }
       
       mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]] <- mantel.partial(bacteriumN_concurrent, bacteriumN2, dist_matrix_contol, method = "pearson", permutations = 999, na.rm = T)
       
       # Store the correlation coefficient and p-value
       cotransmission_BB_cor[bacteriumName2,bacteriumName] <- mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]]$statistic
       cotransmission_BB_pval[bacteriumName2,bacteriumName] <- mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]]$signif
       
     } else {
       mantel_partial_BB_results[[paste0(bacteriumName, '_', bacteriumName2)]] <- 'No concurrent samples'
     }
     
   }
   
 }
 
 ##############################
 # INTERMEDIATE OUTPUT
 ##############################
 # since running the mantel.partial() takes substantial time, raw output should be saved:
 write.table(cotransmission_BB_cor, "03a.RESULTS/Co_transmission_logical_pearson_bacteria_vs_bacteria_corr_raw.txt", sep='\t', quote=F)
 write.table(cotransmission_BB_pval, "03a.RESULTS/Co_transmission_logical_pearson_bacteria_vs_bacteria_pval_raw.txt", sep='\t', quote=F)
 ##############################
 # INTERMEDIATE OUTPUT
 ##############################
 
 colnames(cotransmission_BB_cor) <- check_cotransmission$Host_easy[match(colnames(cotransmission_BB_cor), check_cotransmission$Host_SGB)]
 row.names(cotransmission_BB_cor) <- check_cotransmission$Host_easy[match(row.names(cotransmission_BB_cor), check_cotransmission$Host_SGB)]
  
 colnames(cotransmission_BB_pval) <- check_cotransmission$Host_easy[match(colnames(cotransmission_BB_pval), check_cotransmission$Host_SGB)]
 row.names(cotransmission_BB_pval) <- check_cotransmission$Host_easy[match(row.names(cotransmission_BB_pval), check_cotransmission$Host_SGB)]
  
 cotransmission_BB_cor <- cotransmission_BB_cor[,colSums(is.na(cotransmission_BB_cor)) != nrow(cotransmission_BB_cor)]
 cotransmission_BB_cor <- cotransmission_BB_cor[rowSums(is.na(cotransmission_BB_cor))!=ncol(cotransmission_BB_cor),]

  # since the matrix of correlations is symmetric:   
cotransmission_BB_cor[upper.tri( cotransmission_BB_cor, diag = T)] <- NA
cotransmission_BB_pval[upper.tri( cotransmission_BB_pval,diag = T)] <- NA

# multiple test correction:
cotransmission_BB_flat <- flatten_correlation_matrix(cotransmission_BB_cor, cotransmission_BB_pval)
cotransmission_BB_flat$FDR <- p.adjust(cotransmission_BB_flat$pvalue, method = "BH")

cotransmission_BB_FDR <- cotransmission_BB_pval

 for (i in 1:nrow(cotransmission_BB_flat)) {

   bacteriumName1 <- cotransmission_BB_flat[ i,  'col']
   bacteriumName2 <- cotransmission_BB_flat[ i,  'row']

   cotransmission_BB_FDR[bacteriumName2,bacteriumName1] <- cotransmission_BB_flat[i, 'FDR']

 }
  
cotransmission_BB_cor[is.nan(cotransmission_BB_cor) | !is.finite(cotransmission_BB_cor)] <- NA 
cotransmission_BB_cor[is.nan(cotransmission_BB_cor) | !is.finite(cotransmission_BB_cor)] <- NA
cotransmission_BB_FDR[is.na(cotransmission_BB_cor)] <- NA

indices <- which(cotransmission_BB_cor < -1 | cotransmission_BB_cor > 1, arr.ind = TRUE)

# Print the indices and corresponding values
for (i in 1:nrow(indices)) {
  row <- indices[i, 1]
  col <- indices[i, 2]
  value <- cotransmission_BB_cor[row, col]
  cat("Element at [", row, ",", col, "] is", value, "\n")
}

  # some mantel exact stat artifact I guess:
cotransmission_BB_cor[cotransmission_BB_cor>1] <- 1
cotransmission_BB_cor[cotransmission_BB_cor< (-1)] <- -1
  
pdf('./04.PLOTS/Bacteria_bacteria_co_transmission_corrplot_FDR_0.1_raw.pdf', width=25/2.54, height=22/2.54)
cotransmission_BB_cor[cotransmission_BB_FDR > 0.05] <- 0

 corrplot(cotransmission_BB_cor,
          type = 'lower',
          bg=NA,
          na.label = "X", na.label.col = "grey", tl.col='black')
 
 
 
  dev.off()
##############################
# OUTPUT
##############################
write.table(co_transmission_cor, '03a.RESULTS/Co_transmission.txt', sep='\t', quote=F)
write.table(bgcolors, '03a.RESULTS/Co_transmission_colors.txt', sep='\t', quote=F)

