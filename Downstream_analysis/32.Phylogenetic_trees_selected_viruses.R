setwd('~/Desktop/Projects_2022/NEXT_pilot_FUP/')

#############################################################
# Phylogenetic trees
#############################################################

##############################
# Functions
##############################

##############################
# Loading libraries
##############################
library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)
library(tidytree)
library(tidyverse)
library(MetBrewer)

##############################
# Input data
##############################
VLP_metadata <- read.table('02.CLEAN_DATA/VLP_metadata_final_10_05_2023.txt', sep='\t', header=T)
VLP_metadata$source <- "VLP"
MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"

metadata <- rbind(VLP_metadata[,c("Short_sample_ID","NG_ID","FAM_ID","Timepoint","Type","Individual_ID","Short_sample_ID_bact","Universal_fecal_ID","source")],
                  MGS_metadata[,c("Short_sample_ID","NG_ID","FAM_ID","Timepoint","Type","Individual_ID","Short_sample_ID_bact","Universal_fecal_ID","source")])
metadata$Alex_ID <- paste0(metadata$FAM_ID, '_', metadata$Type, '_', substr(metadata$Short_sample_ID, 1,1), '_', metadata$source, '_', metadata$Timepoint)
row.names(metadata) <- metadata$Alex_ID
metadata$label <- metadata$Alex_ID
metadata$Individual <- substr(metadata$Short_sample_ID, 1,1)
metadata[metadata$Individual=="C",]$Type <- "Infant 1"
metadata[metadata$Individual=="D",]$Type <- "Infant 2"
metadata[metadata$Individual=="A",]$Type <- "Mother"
metadata$tree_ID <- paste0(metadata$FAM_ID, ' ', metadata$Timepoint)

folder = "02.CLEAN_DATA/dist.matrices/"  
file_list = list.files(path=folder, pattern="*.txt")  
virus=lapply(paste0(folder, file_list), function(x) read.table(x, sep='\t', header=T))
names(virus) <- gsub(".dist.txt", '', file_list)

# normalization
virus <- lapply(virus, function(x) {
  x <- x/median( x[upper.tri(x)] )
})

MGS_metadata <- read.table('02.CLEAN_DATA/MGS_metadata_final_10_05_2023.txt', sep='\t', header=T)
MGS_metadata$source <- "MGS"
MGS_metadata$Alex_ID <- paste0(MGS_metadata$FAM_ID, '_', MGS_metadata$Type, '_', substr(MGS_metadata$Short_sample_ID_bact, 1,1), '_', MGS_metadata$source, '_', MGS_metadata$Timepoint)

folder = "01.RAW_DATA/strainphlan_4_distance_matrix/"  
file_list = list.files(path=folder, pattern="*.csv")  
bacterium=lapply(paste0(folder, file_list), function(x) read.csv(x, header=T))
names(bacterium) <- gsub("_dmat_.csv", '', file_list)

# normalization
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


Klimt <- met.brewer('Klimt')
Nizami <- met.brewer('Nizami')

selected_viruses <- read.table('02.CLEAN_DATA/List_viruses_selected_transmission_metadata_with_host.txt', sep='\t', header=T)

##############################
# ANALYSIS
##############################

#### VIRUS L_85266_LS0
tree <- as.phylo(hclust(as.dist(virus[["LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609"]]),method="complete"))

p <- ggtree(tree, options(ignore.negative.edge = T))

p[["data"]]$Type <- metadata$Type[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$FAMILY <- metadata$FAM_ID[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$Sequencing <- metadata$source[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$Timepoint <- metadata$Timepoint[match(p[["data"]]$label, metadata$Alex_ID)]
p[["data"]]$tree_ID <- metadata$tree_ID[match(p[["data"]]$label, metadata$Alex_ID)]

L_85266_LS0_FAMs <- unique(p[["data"]]$FAMILY)
L_85266_LS0_FAMs <- L_85266_LS0_FAMs[!is.na(L_85266_LS0_FAMs)]
L_85266_LS0_label_colors_family <- as.data.frame(L_85266_LS0_FAMs)
L_85266_LS0_label_colors_family$color <- c(Klimt, Nizami[6])

L_85266_LS0_SAMs <- p[["data"]]$tree_ID
L_85266_LS0_SAMs <- L_85266_LS0_SAMs[!is.na(L_85266_LS0_SAMs)]
L_85266_LS0_label_colors <- as.data.frame(L_85266_LS0_SAMs)
L_85266_LS0_label_colors$FAM <- sub("^([^ ]+).*", "\\1", L_85266_LS0_label_colors$L_85266_LS0_SAMs)
L_85266_LS0_label_colors$color <- L_85266_LS0_label_colors_family$color[match(L_85266_LS0_label_colors$FAM, L_85266_LS0_label_colors_family$L_85266_LS0_FAMs)]

pdf('./04.PLOTS/L_85266_LS0_tree.pdf', width=20/2.54, height=18/2.54)
phylo_L_85266_LS0 <- p + geom_tippoint(aes(shape = Type, fill = FAMILY, color=Sequencing), size = 4, stroke=1) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L_85266_LS0_label_colors_family$color, L_85266_LS0_label_colors_family$L_85266_LS0_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.3, size=3) + 
  scale_color_manual(values=setNames(L_85266_LS0_label_colors$color, L_85266_LS0_label_colors$L_85266_LS0_SAMs)) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(aes(shape = Type, fill = FAMILY, colour = Sequencing), size = 5, stroke=1) +
  scale_colour_manual(values = c("VLP"="black", "MGS"="grey")) + 
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  ggtitle("L_85266_LS0") +
  xlim(0, 8)  
dev.off()


#### BACTERIAL HOST (Bacteroides unifromis)
dist_85k_host <- bacterium[["SGB1836_group"]]
# only samples from the same families:
dist_85k_host <- dist_85k_host[grep(paste(L_85266_LS0_FAMs,collapse="|"), 
                                    colnames(dist_85k_host)),
                               grep(paste(L_85266_LS0_FAMs,collapse="|"), 
                                    colnames(dist_85k_host))]
# rarefying the dist matrix so that the phylogenetic tree looks comprehendable:

dist_85k_host <- dist_85k_host[,]

keep_dist_85k_host <- as.data.frame(sort(colnames(dist_85k_host)))
colnames(keep_dist_85k_host) <- "ID"
keep_dist_85k_host[,2:6] <- separate(keep_dist_85k_host, ID, into = c('FAM', 'Type', 'Individual', 'Source', 'Timepoint'), sep = "_")

keep_dist_85k_host$keep <- F

keep_dist_85k_host[gsub("MGS_", '', keep_dist_85k_host$ID) %in%  gsub("VLP_|MGS_", "", colnames(virus[["LN_6A08_VL_306_NODE_3_length_85266_cov_2453.209609"]])), ]$keep <- T



bactree <- as.phylo(hclust(as.dist(dist_85k_host),method="complete"))

bp <- ggtree(bactree, options(ignore.negative.edge = T))

bp[["data"]]$Type <- metadata$Type[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$FAMILY <- metadata$FAM_ID[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$Sequencing <- metadata$source[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$Timepoint <- metadata$Timepoint[match(bp[["data"]]$label, metadata$Alex_ID)]
bp[["data"]]$tree_ID <- metadata$tree_ID[match(bp[["data"]]$label, metadata$Alex_ID)]


Host_L_85266_LS0 <- bp[["data"]]$tree_ID
Host_L_85266_LS0 <- Host_L_85266_LS0[!is.na(Host_L_85266_LS0)]
Host_L_85266_LS0_label_colors_bac <- as.data.frame(Host_L_85266_LS0)
Host_L_85266_LS0_label_colors_bac$FAM <- sub("^([^ ]+).*", "\\1", Host_L_85266_LS0_label_colors_bac$Host_L_85266_LS0)
Host_L_85266_LS0_label_colors_bac$color <- L_85266_LS0_label_colors_family$color[match(Host_L_85266_LS0_label_colors_bac$FAM, L_85266_LS0_label_colors_family$L_85266_LS0_FAMs)]

pdf('./04.PLOTS/Host_L_85266_LS0_tree.pdf', width=15/2.54, height=18/2.54)
phylo_Host_L_85266_LS0 <- bp + geom_tippoint(aes(shape = Type, fill = FAMILY), color = "grey", size = 4, stroke=1) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L_85266_LS0_label_colors_family$color, L_85266_LS0_label_colors_family$L_85266_LS0_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.3, size=3) + 
  scale_color_manual(values=setNames(Host_L_85266_LS0_label_colors_bac$color, Host_L_85266_LS0_label_colors_bac$Host_L_85266_LS0)) +
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  ggtitle("Bacteroides uniformis") + 
  xlim(0, 6) 
dev.off()

#### VIRUS L34922_LS1
tree_L34922_LS1 <- as.phylo(hclust(as.dist(virus[["LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224"]]),method="complete"))

p2 <- ggtree(tree_L34922_LS1, options(ignore.negative.edge = T))

p2[["data"]]$Type <- metadata$Type[match(p2[["data"]]$label, metadata$Alex_ID)]
p2[["data"]]$FAMILY <- metadata$FAM_ID[match(p2[["data"]]$label, metadata$Alex_ID)]
p2[["data"]]$Sequencing <- metadata$source[match(p2[["data"]]$label, metadata$Alex_ID)]
p2[["data"]]$Timepoint <- metadata$Timepoint[match(p2[["data"]]$label, metadata$Alex_ID)]
p2[["data"]]$tree_ID <- metadata$tree_ID[match(p2[["data"]]$label, metadata$Alex_ID)]

L34922_LS1_FAMs <- unique(p2[["data"]]$FAMILY)
L34922_LS1_FAMs <- L34922_LS1_FAMs[!is.na(L34922_LS1_FAMs)]
L34922_LS1_label_colors_family <- as.data.frame(L34922_LS1_FAMs)
L34922_LS1_label_colors_family$color <- c(Klimt, Nizami[6])

L34922_LS1_SAMs <- p2[["data"]]$tree_ID
L34922_LS1_SAMs <- L34922_LS1_SAMs[!is.na(L34922_LS1_SAMs)]
L34922_LS1_label_colors <- as.data.frame(L34922_LS1_SAMs)
L34922_LS1_label_colors$FAM <- sub("^([^ ]+).*", "\\1", L34922_LS1_label_colors$L34922_LS1_SAMs)
L34922_LS1_label_colors$color <- L34922_LS1_label_colors_family$color[match(L34922_LS1_label_colors$FAM, L34922_LS1_label_colors_family$L34922_LS1_FAMs)]

phylo_L34922_LS1 <- p2 + geom_tippoint(aes(shape = Type, fill = FAMILY, color=Sequencing), size = 4, stroke=1) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L34922_LS1_label_colors_family$color, L34922_LS1_label_colors_family$L34922_LS1_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.9, size=3) + 
  scale_color_manual(values=setNames(L34922_LS1_label_colors$color, L34922_LS1_label_colors$L34922_LS1_SAMs)) +
  ggnewscale::new_scale_colour() +
  geom_tippoint(aes(shape = Type, fill = FAMILY, colour = Sequencing), size = 4, stroke=1) +
  scale_colour_manual(values = c("VLP"="black", "MGS"="grey")) + 
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  ggtitle("L34922_LS1") +
  xlim(0, 20)

#### BACTERIAL HOST (Bifidobacterium bifidum)
dist_35k_host <- bacterium[["SGB17256"]]
# only samples from the same families:
dist_35k_host <- dist_35k_host[grep(paste(L34922_LS1_FAMs,collapse="|"), 
                                    colnames(dist_35k_host)),
                               grep(paste(L34922_LS1_FAMs,collapse="|"), 
                                    colnames(dist_35k_host))]

# only concurrent samples from the same families:
concurrent <- grep(paste(gsub('VLP','MGS',colnames(virus[["LN_4F02_VL_264_NODE_17_length_34922_cov_30951.540224"]])),collapse="|"), 
     colnames(dist_35k_host), value = T)

FAM0234_bifidum <- grep('FAM0234', colnames(dist_35k_host), value=T)
FAM0432_bifidum <- grep('FAM0432', colnames(dist_35k_host), value=T)

dist_35k_host_samples <- unique(c(concurrent,FAM0234_bifidum, FAM0432_bifidum))

dist_35k_host <- dist_35k_host[dist_35k_host_samples,
                               dist_35k_host_samples]

bactree2 <- as.phylo(hclust(as.dist(dist_35k_host),method="complete"))

bp2 <- ggtree(bactree2, options(ignore.negative.edge = T))

bp2[["data"]]$Type <- metadata$Type[match(bp2[["data"]]$label, metadata$Alex_ID)]
bp2[["data"]]$FAMILY <- metadata$FAM_ID[match(bp2[["data"]]$label, metadata$Alex_ID)]
bp2[["data"]]$Sequencing <- metadata$source[match(bp2[["data"]]$label, metadata$Alex_ID)]
bp2[["data"]]$Timepoint <- metadata$Timepoint[match(bp2[["data"]]$label, metadata$Alex_ID)]
bp2[["data"]]$tree_ID <- metadata$tree_ID[match(bp2[["data"]]$label, metadata$Alex_ID)]

Host_L34922_LS1 <- bp2[["data"]]$tree_ID
Host_L34922_LS1 <- Host_L34922_LS1[!is.na(Host_L34922_LS1)]
Host_L34922_LS1_label_colors_bac <- as.data.frame(Host_L34922_LS1)
Host_L34922_LS1_label_colors_bac$FAM <- sub("^([^ ]+).*", "\\1", Host_L34922_LS1_label_colors_bac$Host_L34922_LS1)
Host_L34922_LS1_label_colors_bac$color <- L34922_LS1_label_colors_family$color[match(Host_L34922_LS1_label_colors_bac$FAM, L34922_LS1_label_colors_family$L34922_LS1_FAMs)]


phylo_Host_L34922_LS1 <- bp2 + geom_tippoint(aes(shape = Type, fill = FAMILY), color = "grey", size = 4, stroke=1) +
  scale_shape_manual(values=c("Mother" = 21, "Infant 1"=24, "Infant 2"=23)) + 
  scale_fill_manual(values=setNames(L34922_LS1_label_colors_family$color, L34922_LS1_label_colors_family$L_85266_LS0_FAMs)) + 
  guides(color="none") +
  geom_tiplab(aes(label=tree_ID, color=tree_ID), offset=0.02, size=3) + 
  scale_color_manual(values=setNames(Host_L34922_LS1_label_colors_bac$color, Host_L34922_LS1_label_colors_bac$Host_L34922_LS1)) +
  geom_treescale() +
  guides(shape = guide_legend(override.aes = list(color = NULL)),
         fill = guide_legend(override.aes = list(shape = 21, color = "black"))) +
  ggtitle("Bifidobacterium bifidum") + 
  xlim(0, 0.8) 


##############################
# OUTPUT
##############################
###### FOR VISUALIZATION #####
saveRDS(p, file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_L_85266_LS0.rds")
saveRDS(bp, file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_Host_L_85266_LS0.rds")
write.table(L_85266_LS0_label_colors, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_L_85266_LS0_label_colors.txt", sep='\t', row.names=F)
write.table(L_85266_LS0_label_colors_family, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_L_85266_LS0_family_colors.txt", sep='\t', row.names = F)
write.table(Host_L_85266_LS0_label_colors_bac, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5B_Host_L_85266_LS0_label_colors.txt", sep='\t', row.names=F)

saveRDS(p2, file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_L34922_LS1.rds")
saveRDS(bp2, file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_Host_L34922_LS1.rds")
saveRDS(bp2, file = "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_Host_L34922_LS1_concurrent.rds")
write.table(L34922_LS1_label_colors, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_L34922_LS1_label_colors.txt", sep='\t', row.names=F)
write.table(L34922_LS1_label_colors_family, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_L34922_LS1_label_family_colors.txt", sep='\t', row.names = F)
write.table(Host_L34922_LS1_label_colors_bac, "02.CLEAN_DATA/PREPARED_DATA_FOR_PLOTS/Fig5C_Host_L34922_LS1_label_colors.txt", sep='\t', row.names=F)

