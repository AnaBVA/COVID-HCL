library(rhdf5)
library(SingleCellExperiment)
library(HDF5Array)

file <- "data/17727365"
#hcl <- h5read(file,"X")

## Creating a SingleCellExperiment
# Open h5File
h5 <- HDF5Array(file,"X")
dim(h5)

h <- H5Fopen(file)

sce <- SingleCellExperiment(
    list(counts = h5),
    rowData = h$var,
    colData = h$obs
)

h5closeAll()

# Open meta data
meta <- read.csv("data/HCL_Fig1_cell_Info.csv",header = T)

cellInfo <- data.frame(meta)
cellTypeColumn <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))

saveRDS(cellInfo, file="int/cellInfo.Rds")



### Motif database
# featherURL <- "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather"
# download.file(featherURL, destfile=basename(featherURL))

motifRankings <- importRankings("hg19-tss-centered-10kb-7species.mc9nr.feather")
data(motifAnnotations_hgnc)

org="hgnc" # or hgnc, or dmel
dbDir="cisTarget_databases" # RcisTarget databases location
myDatasetTitle="SCENIC in HCL data" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)





