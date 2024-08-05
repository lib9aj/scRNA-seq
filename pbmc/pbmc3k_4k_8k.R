
library(Seurat)

pbmc3k.data <- Read10X (data.dir = "/Users/lib9aj/git/scRNA-seq/pbmc/pbmc3k/filtered_gene_bc_matrices/hg19")
pbmc4k.data <- Read10X (data.dir = "/Users/lib9aj/git/scRNA-seq/pbmc/pbmc4k/filtered_gene_bc_matrices/GRCh38")
pbmc8k.data <- Read10X (data.dir = "/Users/lib9aj/git/scRNA-seq/pbmc/pbmc8k/filtered_gene_bc_matrices /GRCh38")

pbmc3k <- CreateSeuratObject(counts = pbmc3k.data, project = "PBMC3K")
pbmc4k <- CreateSeuratObject(counts = pbmc4k.data, project = "PBMC4K")
pbmc8k <- CreateSeuratObject(counts = pbmc8k.data, project = "PBMC8K")

pbmc3k
pbmc4k
pbmc8k

###Merging Two Seurat Objects

pbmc.combined <- merge (pbmc4k, y = pbmc8k, add.cell.ids = c ("4K", "8K"), project = "PBMC12K")
pbmc.combined

# notice the cell names now have an added identifier
head(colnames(pbmc.combined))

table (pbmc.combined $orig.ident)
table (pbmc.combined $nCount_RNA)


###Merging More Than Two Seurat Objects
#using the 4K and 8K PBMC datasets as well as our previously computed Seurat object from the 2,700 PBMC tutorial 
#(loaded via the SeuratData package）
library(SeuratData)
InstallData("pbmc3k")
pbmc3k <- LoadData("pbmc3k", type = "pbmc3k.final")
pbmc3k
pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = c("3K", "4K", "8K"), project = "PBMC15K")
pbmc.big

head(colnames(pbmc.big))
tail(colnames(pbmc.big))

unique (sapply (X = strsplit(colnames(pbmc.big), split = "_"), FUN = "[", 1))

table(pbmc.big$orig.ident)


