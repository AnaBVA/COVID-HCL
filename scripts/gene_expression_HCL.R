############################################################
############################################################
# ACE2 and TPRMSS2 gene expression in single cell data
############################################################
## AIM:
############################################################
#
# This is script is to plot ACE2 and TPRMSS2 gene expression
# from sc data using the Human Cell Landscate data.
#
############################################################
## Downloaded data from The Human Cell Landscape (HCL)
############################################################
#
# Paper: “Construction of a human cell landscape at single-cell level” (Han et al. 2020)
# Web interface: https://db.cngb.org/HCL/landscape.html
# Data: https://figshare.com/articles/HCL_DGE_Data/7235471
# Github: https://github.com/ggjlab/HCL
# R packge: scHCL  https://github.com/ggjlab/scHCL/
#
######################### Note
#
# For the analysis we used:
# 1) “HCL_Fig1_adata.h5ad” file
# 2) "HCL_Fig1_cell_Info.xls" file
#
############################################################
## Setup
############################################################
#
# This script was run in the "DNA" cluster from LAVIS at UNAM
#  ssh -X aaltamirano@dna.lavis.unam.mx
# qrsh -l h_vmem=130G
# module load r/3.6.1
# /home/amedina/amedinalab/aaltamirano/scRNAseq/COVID19
#
############################################################
# Libraries
############################################################

library(rhdf5)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(igraph)
library(stringr)

############################################################
# Data
############################################################
# Meta data
meta <- read.csv("data/HCL_Fig1_cell_Info.csv",header = T)

# H5 file location
file <- "data/17727365"

# file structure
h5ls(file)

# Open file into R env
h <- H5Fopen("data/17727365")

# ACE2 and TMPRSS2 gene possition
p.ace2 <- which(h$var$index == "ACE2")
p.tmprss2 <- which(h$var$index == "TMPRSS2")

# Index, batch, tissue, number of genes and number of counts per cell in a df
obs <- h$obs

# Save values
obs$ace2 <- h$X[p.ace2,]
obs$tmprss2 <- h$X[p.tmprss2,]

# Trim index name
obs$cellnames <- gsub("\\-[0-9].*","",obs$index)

# Merge obs with expression values and meta data info
all <- merge(obs,meta,by="cellnames")

h5closeAll()

#write.csv(all,"all.txt",row.names = F, quote = F)
#all <- read.csv("all.txt",header = T)

############################################################
# ACE2
############################################################

df <-all %>%
    filter(ace2 != 0) %>%
    arrange(desc(ace2, FUN = median)) %>%
    group_by(celltype) %>%
    #filter(sample=="AdultLung"|str_detect(sample,"Fetal")) %>%
    filter(max(ace2)>=2)


p <- ggplot(df,aes(x=factor(celltype, levels = unique(celltype)),y =ace2)) +
        geom_boxplot() +
        #geom_jitter(data=. %>% filter(sample =="AdultArtery"),aes(colour=sample, alpha=0.5)) + # data=.%>% filter(sample =="AdultLung")
        geom_jitter(aes(colour=sample, alpha=0.5))+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              text = element_text(size = 30))+
        #scale_color_manual(values = colorRampPalette(brewer.pal(12, "Dark2"))(length(unique(all$sample)))) +
        labs(title= "ACE2 expression", y="Counts in ln(CPM/100 + 1) ", x = "Celltype")


pdf(str_c("fig/ACE2_celltypes_",Sys.Date(),".pdf"),width = 26,height = 15)
p
dev.off()

# pp <- ggplotly(p)
# saveRDS(pp, "ACE2.rds")


# n <-all %>%
#     group_by(sample) %>%
#     filter(sample=="AdultLung") %>%
#     filter(ace2 != 0)
#
# n$cell <- gsub('[()]',"",n$cell)
#
# network <- graph_from_data_frame(d=n, directed=F)
# network <- graph_from_literal(strsplit(n$cells[1],"-"))
# network <- graph_from_literal("AT2 cell":"Macrophage":"Endothelial cell (endothelial to mesenchymal transition)":"M2 Macrophage":"Smooth muscle cell":"Basal cell":"T cell":"B cell (Plasmocyte)":"Stratified epithelial cell":"Monocyte":"Endothelial cell (APC)":"Dendritic cell":"Mast cell":"Fibroblast":"Epithelial cell (intermediated)":"Endothelial cell":"Sinusoidal endothelial cell":"Chondrocyte":"Stromal cell":"CB CD34+":"Fetal chondrocyte":"Fetal epithelial progenitor":"Antigen presenting cell (RPS high)":"Neutrophil (RPS high)":"B cell":"Goblet cell":"Fetal mesenchymal progenitor":"Fetal fibroblast":"Fetal stromal cell":"Epithelial cell":"Kidney intercalated cell":"Loop of Henle":"Gastric endocrine cell":"Myeloid cell":"Mesothelial cell":"Enterocyte progenitor":"Fasciculata cell")
#
# deg <- n$ace2 *10
#
# pdf("ACE2_network.pdf")
# plot(network,vertex.size=deg)
# dev.off()

# a <- all %>%
#     group_by(cluster) %>%
#     summarise(mean=mean(ace2),
#               mean_without_0 = mean(which(ace2!=0)),
#               ncells_without_0 = sum(ace2!=0),
#               ncells=n(),
#               tissue =paste(sample, collapse = ","),
#               celltype=paste(celltype, collapse = ","),
#               donor = paste(donor, collapse = ",")
#               ) %>%
#     arrange(desc(mean_without_0))
#
# write.csv(a,"ace2_hcl.txt",row.names = F, quote = F)

############################################################
# TMPRSS2
############################################################

df <-all %>%
    filter(tmprss2 != 0) %>%
    arrange(desc(tmprss2, FUN = median)) %>%
    group_by(celltype) %>%
    filter(max(tmprss2)>=2)
    #mutate(tissue = paste(unique(sample), collapse="-"))


p <- ggplot(df,aes(x=factor(celltype, levels = unique(celltype)),y =tmprss2)) +
    geom_boxplot() +
    #geom_jitter(data=. %>% filter(sample =="AdultArtery"),aes(colour=sample, alpha=0.5)) + # data=.%>% filter(sample =="AdultLung")
    geom_jitter(aes(colour=sample, alpha=0.5))+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size = 30))+
    labs(title= "Tmprss2 expression", y="Counts in ln(CPM/100 + 1) ", x = "Celltype")


pdf(str_c("fig/TMPRSS2_celltypes_",Sys.Date(),".pdf"),width = 26,height = 15)
p
dev.off()


# t <- obs %>%
#     group_by(tissue) %>%
#     summarise(mean=mean(tmprss2),
#               mean_without_0 = mean(which(tmprss2!=0)),
#               ncells_without_0 = sum(tmprss2!=0),
#               ncells=n(),
#               ntissue =first(gsub("\\_.*","",index))
#     ) %>%
#     arrange(desc(mean_without_0))
#
#
# write.csv(t,"tmprss2_hcl.txt",row.names = F, quote = F)

