#setwd("../TissueTimer/DaphneCode")
library(data.table)

####################
#Step 1: read in Christoph expression data
load.data <- function() {
    
    c.db <<- fread("originalExpressionData/christoff_data.txt")
    c.meta <<- fread("originalExpressionData/christoff_meta.csv")
    
}

load.data()

db <- c.db
meta <- c.meta
colnames(meta)=paste(meta[1,])
meta=meta[-1,]

####################
#Step 2: make files for cibersort

load.data()
#setwd("../../../DaphneCode")
db <- c.db #reloading data because CIBERTSORT wants linear, non-log transformed TPM
data <- db[, names(db) %in% c("gene", meta$sample), with = FALSE]

keepers <- c("root", "hypocotyl", "stem", "cotyledons", "leaf", "apex")  #The only tissues we care about

meta <- meta[clean.tissue %in% keepers]

#selecting roots
meta <- meta[clean.tissue != "root" | subtissue == "root_whole"]

#proper label for hypocotyl as independent category
#done by default

#selecting stems 
meta[clean.tissue == 'stem', subtissue := "stem"]

#selecting leaves
meta <- meta[clean.tissue != "leaf" | !grepl("petiole", subtissue)]
meta <- meta[clean.tissue != "leaf" | !grepl("leaf tip", subtissue)]
meta[clean.tissue == "cotyledons", clean.tissue := 'leaf'] #merging cotyledons with leaves
meta[clean.tissue == "leaf", subtissue := "leaf"]

#selecting apex
meta <- meta[clean.tissue != "apex" | !grepl("clv", subtissue)]

#rewrangling flowers - only selecting those that are whole flowers 
#meta <- meta[clean.tissue != "flower" | grepl('[0-9]+$', subtissue)]
#meta[clean.tissue == "flower", subtissue := "whole"]

setorder(meta, clean.tissue, subtissue, time)
meta[,tissue := paste(clean.tissue, subtissue, time, sep = "_")]

tmp.data <- data[,meta$sample, with = F]
data <- cbind(gene = data[,gene], tmp.data)

###Save table that will be used to find markers using cibersort
write.table(data, "data_to_cibersort.txt", quote = F, row.names = F, col.names = T, sep = "\t")

####################
# Step 3: Process phenotype data
pheno <- data.frame(dcast(tissue ~ sample, data = meta, fill = 2))
rownames(pheno) <- pheno$tissue
pheno$tissue <- NULL
pheno <- ifelse(pheno == 2, 2, 1)
out.pheno <- pheno[, meta$sample]

all(colnames(out.pheno) == names(data)[2:ncol(data)]) # must be true

write.table(out.pheno, "classes_to_cibersort.txt", quote = F, row.names = T, col.names = F, sep = "\t")

## the above line generates the comparison for all classes to all classes. now we'll consider the 
## comparison where within-tissue comparisons are ignored
## put 0s in every spot within a tissue of the same type

pheno <- data.frame(dcast(tissue ~ sample, data = meta, fill = 2))
rownames(pheno) <- pheno$tissue
pheno$tissue <- NULL
pheno <- ifelse(pheno == 2, 2, 1)
out.pheno <- pheno[, meta$sample]

for (i in 1:nrow(out.pheno)) {
    for (j in 1:ncol(out.pheno)) {
        if (out.pheno[i,j] == 2) {
            if (meta[sample == colnames(out.pheno)[j], clean.tissue] == strsplit(rownames(out.pheno)[i], "_")[[1]][1]){
                out.pheno[i,j] <- 0
            }
        }
    }
}

all(colnames(out.pheno) == names(data)[2:ncol(data)]) # must be true

write.table(out.pheno, "classes_to_cibersort_ignore_in_tissue.txt", quote = F, row.names = T, col.names = F, sep = "\t")



## now we'll ignore within-timepoint comparisons 

pheno <- data.frame(dcast(tissue ~ sample, data = meta, fill = 2, value.var = "tissue"))
rownames(pheno) <- pheno$tissue
pheno$tissue <- NULL
pheno <- ifelse(pheno == 2, 2, 1)
out.pheno <- pheno[, meta$sample]

meta[, mod.time := gsub("+", "", time, fixed = TRUE)]

for (i in 1:nrow(out.pheno)) {
    for (j in 1:ncol(out.pheno)) {
        if (out.pheno[i,j] == 2) {
            if (meta[sample == colnames(out.pheno)[j], time] == rev(strsplit(rownames(out.pheno)[i], "_")[[1]])[1]){
                out.pheno[i,j] <- 0
            }
        }
    }
}


all(colnames(out.pheno) == names(data)[2:ncol(data)]) # must be true

write.table(out.pheno, "classes_to_cibersort_ignore_in_time.txt", quote = F, row.names = T, col.names = F, sep = "\t")

#renaming to remove second _ 

meta[, inter.tissue := paste(clean.tissue, subtissue, sep = "_")]
meta[, renamed := paste(inter.tissue, gsub('+', '', gsub("d", '',time), fixed = TRUE), sep = "DATE")]


