require(ggplot2)
#require(LaCroixColoR)
require(data.table)
#require(lettercase)
## testing code
source("tissue_timer.R")
#example boundary list for some of our data
boundary.list <- list(apex_inflorescence = c(14, 28), apex_vegetative = c(NA,21), flower_whole = c(14, NA), 
                      hypocotyl_hypocotyl = c(NA, 15), leaf_leaf = c(NA,NA), root_root_whole = c(NA,NA), stem_stem = c(NA, NA))



#mp <- "~/Documents/Alan Turing/Projects/data/signatures_from_cibersort.txt"
#tp <- "~/Documents/Alan Turing/Projects/data/christoff_data.txt"
mp <- "signatures/scaled_time_signatures.txt"
tp <- "originalExpressionData/christoff_data.txt"

## testing on training data
outputs <- main(marker.path = mp, test.samples = tp, bd.list = boundary.list, start.in = 0, stop.in = 28, error.fun = 'mono')

# takes result of main function
output <- outputs[, .SD[which.max(percent), .(tissue, time)], by = sample.name]
setnames(output, c('sample', 'predicted.tissue', 'predicted.time'))

merged.test <- output[meta, on = "sample"]

merged.test[clean.tissue == "apex", 
            clean.tissue := paste(clean.tissue, ifelse(grepl("inf", subtissue), "inflorescence", "vegetative"), sep = "_")]

mapped_tissue_dict <- unique(merged.test$predicted.tissue)
true_tissue_dict <- unique(merged.test$clean.tissue)

names(mapped_tissue_dict) <- c("root", 'hypocotyl', 'stem', 'leaf', 'apex_vegetative', "apex_inflorescence", 'flower')
mapping <- data.table(pred.names <- mapped_tissue_dict, corr.names <- names(mapped_tissue_dict))
setnames(mapping, c("predicted.tissue", 'named.predicted'))

merged.test <- merged.test[mapping, on = 'predicted.tissue']
merged.test[, real.time := gsub("d*\\+*", "", time)]

merged.test[, clean.tissue := factor(clean.tissue, levels = c('root', 'hypocotyl', 'stem', 'apex_vegetative', 'leaf', 'cotyledons', 'flower', 'apex_inflorescence', 'silique', 'seed'))]
merged.test[, named.predicted := factor(named.predicted, levels = c('root', 'hypocotyl', 'stem', 'apex_vegetative', 'leaf', 'cotyledons', 'flower', 'apex_inflorescence', 'silique', 'seed'))]
breakdown <- merged.test[, c('sample', 'clean.tissue', 'real.time', 'named.predicted', 'predicted.time'), with = F]

###We don't care about flower, silique, or seed
breakdown_sub=breakdown[which(as.character(breakdown[,'clean.tissue'][[1]]) %in% c("root", "hypocotyl", "stem", "apex_vegetative", "leaf", "cotyledons", "apex_inflorescence")),]

####How many correct predictions?
#correct
table(as.character(breakdown_sub[,'clean.tissue'][[1]])[which(breakdown_sub[,'clean.tissue']==breakdown_sub[,"named.predicted"])])
#incorrect
table(as.character(breakdown_sub[,'clean.tissue'][[1]])[which(breakdown_sub[,'clean.tissue']!=breakdown_sub[,"named.predicted"])])

###what were the 6 mistakes?
breakdown_sub[,c('clean.tissue', 'named.predicted', 'real.time')][which(breakdown_sub[,'clean.tissue']!=breakdown_sub[,"named.predicted"])]

#####Now, let's look at the age differences
breakdown_sub[1:72,c("real.time", "predicted.time")]
offBy=as.numeric(breakdown_sub[,c("real.time")][[1]])-as.numeric(breakdown_sub[,c("predicted.time")][[1]])


relevantTissues=c("apex_vegetative", "apex_inflorescence","cotyledons", "hypocotyl", "leaf", "root", "stem")
relevantTissueColors=c("darkmagenta", "lightpink", "lightgreen", "grey", "forestgreen", "orange", "cornflowerblue")
names(relevantTissueColors)=relevantTissues


pdf(paste("Figure_S2_TissueTimer_predictionsPureTissue_scatter.pdf", sep=""), width=3, height=3,pointsize = 10) 
par(mfrow=c(1,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.5)

set.seed(1)
jittered_x=jitter(as.numeric(breakdown_sub[,c("real.time")][[1]]), factor=2)
jittered_y=jitter(as.numeric(breakdown_sub[,c("predicted.time")][[1]]), factor=2)
plot(jittered_x,jittered_y,
     xlab="real age (days) -jittered", ylab="predicted age (days) -jittered", pch=19, col=relevantTissueColors[as.character(breakdown_sub[,c("clean.tissue")][[1]])])
abline(c(0,1))

legend(27,20, names(relevantTissueColors), col=relevantTissueColors, pch=19)

dev.off()

pdf(paste("Figure_2a_TissueTimer_predictionsPureTissue.pdf", sep=""), width=3, height=4,pointsize = 10) 
par(mfrow=c(3,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.5)

diff=as.numeric(breakdown_sub[,c("real.time")][[1]])-as.numeric(breakdown_sub[,c("predicted.time")][[1]])
tissues=as.character(breakdown_sub[,c("clean.tissue")][[1]])
plot(table(diff[which(tissues=="leaf")]), col="forestgreen", xlab="real age - predicted age", ylab="number of samples", xlim=c(-20, 20), main="leaf, sd=3.65")
plot(table(diff[which(tissues=="root")]), col="orange", xlab="real age - predicted age", ylab="number of samples", xlim=c(-20, 20), main="root, sd=2.03")
plot(table(diff[which(tissues!="root" & tissues!="leaf")]), col="darkgrey", xlab="real age - predicted age", ylab="number of samples", xlim=c(-20, 20), main="other tissue, sd=7.90")

dev.off()

sd(diff[which(tissues=="leaf")])
sd(diff[which(tissues=="root")])
sd(diff[which(tissues!="root" & tissues!="leaf")])





#####Repeat this analysis for trava database

meta_trava=read.table("originalExpressionData/metaTravaDB.csv", header=T, stringsAsFactors=F, sep=",")
rownames(meta_trava)=meta_trava[,2]
db_trava = read.table("originalExpressionData/TPM.travadb.txt", header=T, sep="\t", stringsAsFactors = F)

meta_trava[,"clean.tissue"]=tolower(meta_trava[,'clean.tissue'])
colnames(meta_trava)[2]="sample"
trava.keep <- c("cotyledon", 'hypocotyl','internode','leaf','root')

meta_trava <- meta_trava[meta_trava[,"clean.tissue"] %in% trava.keep, ]
setnames(db_trava, 'gene_id', 'gene')

db_trava <- db_trava[, which(names(db_trava) %in% c('gene', rownames(meta_trava)))] 

fwrite(db_trava, "trava_test.txt")
tp <- "trava_test.txt"
outputs_trava <- main(marker.path = mp, test.samples = tp, bd.list = boundary.list, start = 0, stop = 28)
output_trava <- outputs_trava[, .SD[which.max(percent), .(tissue, time)], by = sample.name]
setnames(output_trava, c('sample', 'predicted.tissue', 'predicted.time'))

merged.test_trava <- output_trava[meta_trava, on = "sample"]
merged.test_trava[grepl('inflorescence', subtissue), clean.tissue := "apex_inflorescence"]
merged.test_trava[grepl('internode', clean.time), clean.tissue := "stem"]
replace <- c('hypocotyl', 'apex_vegetative', 'apex_inflorescence', 'leaf', 'flower', 'stem', 'root')
names(replace) <- c("hypocotyl_hypocotyl","apex_vegetative", "apex_inflorescence",  "leaf_leaf", "flower_whole", "stem_stem", "root_root_whole")

merged.test_trava[, predicted.tissue := replace[predicted.tissue]]

#in clean.tissue there is 'internode' which is a kind of stem
merged.test_trava[which(as.character(merged.test_trava[,"clean.tissue"][[1]])=="internode"),"clean.tissue"]="stem"

#compare tissue accuracy
#correct
table(as.character(merged.test_trava[,'clean.tissue'][[1]])[which(merged.test_trava[,'clean.tissue']==merged.test_trava[,"predicted.tissue"])])
#incorrect
table(as.character(merged.test_trava[,'clean.tissue'][[1]])[which(merged.test_trava[,'clean.tissue']!=merged.test_trava[,"predicted.tissue"])])

###what were the mistakes
merged.test_trava[which(merged.test_trava[,'clean.tissue']!=merged.test_trava[,"predicted.tissue"]),]

#two cotyledons mislabelled as 6 and 8 day old leaves
#two hypocotyls mislabelled as 4-day old leaves
#two stems (senescent internode) mislabelled as leaves 
#two leaf petioles of young leaves mislabelled as vegetative apex

pdf(paste("Figure_S2_TissueTimer_predictionsTravaDB.pdf", sep=""), width=3, height=4,pointsize = 10) 
par(mfrow=c(2,1));
par(mar=c(4, 4, 1, 1)+0.1);
par(cex=0.5)

#what are the ages of leaves?
merged.test_trava[which(merged.test_trava[,'clean.tissue']=="leaf"),]

agesPred_leaf=lapply(c("Young", "Intermediate", "Mature", "Scenescent"), function(i){
    as.numeric(merged.test_trava[which(merged.test_trava[,'time']==i),"predicted.time"][[1]])
})
names(agesPred_leaf)=c("Young (n=4)", "Intermediate (n=10)", "Mature (n=8)", "Scenescent (n=6)")
boxplot(agesPred_leaf, xlab="leaf maturity", ylab="predicted age (days)")

merged.test_trava[which(merged.test_trava[,'clean.tissue']=="root")]
agesPred_root=lapply(c("seedling.root", "root.apex", "root"), function(i){
    as.numeric(merged.test_trava[which(merged.test_trava[,'suggested.tissue']==i),"predicted.time"][[1]])
})
names(agesPred_root)=c("root seeding, 1d (n=2)", "root apex, 7d (n=2)", "root no apex, 7d (n=2)")
boxplot(agesPred_root, xlab="root sample", ylab="predicted age (days)")

dev.off()




