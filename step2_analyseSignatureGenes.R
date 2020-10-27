library(gplots)
library(pheatmap)
sig_tissues=read.table("signatures/scaled_time_signatures_ignore_time.txt", sep=",", head=T, row.names=1)
sig_time=read.table("signatures/scaled_time_signatures_ignore_tissue.txt", sep=",", head=T, row.names=1)
sig_both=read.table("signatures/scaled_time_signatures.txt", sep=",", head=T, row.names=1)

#Fig 1A: Venn diagram
genes_tissues=paste(rownames(sig_tissues), "")
genes_time=paste(rownames(sig_time), "")
genes_both=paste(rownames(sig_both), "")

length(genes_tissues)
length(genes_time)
length(genes_both)

#middle
length(which(genes_tissues %in% genes_time & genes_tissues %in% genes_both))

#pairwise intersections
length(which(genes_tissues %in% genes_time & !(genes_tissues %in% genes_both)))
length(which(genes_tissues %in% genes_both & !(genes_tissues %in% genes_time)))
length(which(genes_time %in% genes_both & !(genes_time %in% genes_tissues)))

#unique genes
length(which(!(genes_tissues %in% genes_time) & !(genes_tissues %in% genes_both)))
length(which(!(genes_both %in% genes_tissues) & !(genes_both %in% genes_time)))
length(which(!(genes_time %in% genes_both) & !(genes_time %in% genes_tissues)))


genes_tissues[which(!(genes_tissues %in% genes_time) & !(genes_tissues %in% genes_both))]
#AT2G38540: LIPID TRANSFER PROTEIN 1 (most expressed in mature flower according to eFP)
#AT5G20630: GERMIN 3 (most expressed in leaf and young flower)
#AT5G54370: Late embryogenesis abundant (most expressed in root)

genes_time[which(!(genes_time %in% genes_both) & !(genes_time %in% genes_tissues))]
#Highlughts include:
#AT2G26330: Erecta: Involved in specification of organs originating from the shoot apical meristem
#several members of light harvesting complex LHCA1,2,3,4,4.2 PSI-P, PSAD-2, PSAH, PSBR
#AT4G03210: XYLOGLUCAN ENDOTRANSGLUCOSYLASE/HYDROLASE: encodes a member of xyloglucan endotransglucosylase/hydrolases (XTHs) that catalyze the cleavage and molecular grafting of xyloglucan chains function in loosening and rearrangement of the cell wall.




###from gProfiler-- prepare files

write.table(genes_both, file="genes_in_sigMatrix.txt", row.names = F, col.names = F, quote=F)

###make heatmaps and/or pca plots

heatmap.2(-as.matrix(sig_both), trace='n', scale='n')
#use original cristoff table
load.data <- function() {
    c.db <<- read.table("originalExpressionData/christoff_data.txt", header=T, row.names=1, stringsAsFactors = F)
    c.meta <<- fread("originalExpressionData/christoff_meta.csv")
}

load.data()

db <- c.db
meta <- c.meta
colnames(meta)=paste(meta[1,])
meta=meta[-1,]

relevantTissues=c("apex", "cotyledons", "hypocotyl", "leaf", "root", "stem")
relevantTissueColors=c("darkmagenta", "lightgreen", "yellow", "forestgreen", "orange", "cornflowerblue")

relevantAges=c("5d", "7d", "10d", "14d", "17d", "21d", "21d+", "24d", "25d", "28d", "28d+", "35d")
relevantAgeColors=grey.colors(12)
names(relevantAgeColors)=relevantAges

names(relevantTissueColors)=relevantTissues
db_sub=db[rownames(sig_both), which(as.character(meta[,"clean.tissue"][[1]])  %in% relevantTissues)]

ann_colors= list("clean.tissue"=relevantTissueColors, "time"=relevantAgeColors)

annots=data.frame(meta[,c("clean.tissue", "time")])
rownames(annots)=as.character(meta[,"sample"][[1]])

pheatmap(log(db_sub+1), legend=T, annotation_col=annots, annotation_colors=ann_colors, labels_col="sample", labels_row="marker gene", file="Figure_1_TissueTimer_heatmapOfMarkers.pdf",cutree_cols=3, show_rownames=F, show_colnames=F)

#sub-plot: only include markers that were unique to time, not tissue
temp=trimws(genes_time[which(!(genes_time %in% genes_both) & !(genes_time %in% genes_tissues))])
db_sub2=db[temp, which(as.character(meta[,"clean.tissue"][[1]])  %in% relevantTissues)]
db_sub2=log(db_sub2+1)
pheatmap(db_sub2, legend=T, annotation_col=annots, annotation_colors=ann_colors,cutree_cols=3, show_rownames=T, show_colnames=F)



