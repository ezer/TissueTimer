####preliminary stuff to set up, so this works as an independent peice of code:
library(MASS)
props=seq(0, 1, 0.1)

propsByTissue=unlist(sapply(props, function(leaf){
    tot=leaf
    sapply(props[which(props<(1-tot))], function(stem){
        tot=tot+stem
        lapply(props[which(props<(1-tot))], function(root){
            c(leaf, stem, root, 1-root-stem-leaf)   
        })
    })
}))

propsByTissue=t(matrix(propsByTissue, nrow =4))  

propsByTissue_Sample=paste("Sample", apply(propsByTissue, 1, function(i){ paste(i, collapse="_")}), sep="_")



colnames(propsByTissue)=c("leaf", "stem", "root", "apex")

#for every tissue ratio being investigated, generatea a mean tissue value:
mp <- "signatures/scaled_time_signatures.txt"
markers <- load.markers(mp) #gets markers that we want. these are currently generated via CIBERSORT

#make fake 'marker' files, where the gene expression values are tpms from single samples, instead of averages
cibersort_input=read.table(file="cibersortInputs/data_to_cibersort.txt")
cibersort_labels=read.table(file="cibersortInputs/classes_to_cibersort.txt")





#######  We have three strategies of modelling noise:

bd.list <- list(apex_inflorescence = c(14, 28), apex_vegetative = c(NA,21), flower_whole = c(14, NA), 
                hypocotyl_hypocotyl = c(NA, 15), leaf_leaf = c(NA,NA), root_root_whole = c(NA,NA), stem_stem = c(NA, NA))
start.in=2
stop.in=14
step.in=2

####Method 1. Use lots of different randomly sampled permutations for linear combinations
set.seed(1)
fakeMarkerSets=lapply(c(1:10), function(i){
    a=apply(cibersort_labels, 1, function(j){
        ids=sample(which(as.character(j)=="1"), 1)
        cibersort_input[,ids]
    })  
    a=cbind(as.character(cibersort_input[,1]), a)
    colnames(a)=c("gene", cibersort_labels[,1])
    a
})
#add blank flowers:
fakeMarkerSets=lapply(fakeMarkerSets, function(i){
    temp=cbind(i[,1:6], rep(0, length(i[,1])), i[,7:18])
    colnames(temp)=colnames(markers)
    temp
})
#add blank flowers:
fakeMarkerSets=lapply(fakeMarkerSets, function(i){
    temp=cbind(i[,1:6], rep(0, length(i[,1])), i[,7:18])
    colnames(temp)=colnames(markers)
    temp
})
#inter.marks <- interpolate.markers(markers, boundary.list = bd.list, start = start.in, stop = stop.in) #interpolates onto the timescale we're intersted in in
markerGenes=as.character(markers[,1][[1]])
set.seed(1)

colnames(propsByTissue)=c("leaf_leaf", "stem_stem", "root_root_whole", "apex_vegetative")
sim1 <- apply(propsByTissue, 1,  function(j){
    lapply(c(1:10), function(i){
        temp=fakeMarkerSets[[i]]
        ids=which(as.character(fakeMarkerSets[[i]][,1]) %in% markerGenes)
        temp=temp[ids,]
        temp=interpolate.markers(data.table(temp), boundary.list = bd.list, start = start.in, stop = stop.in) #interpolates onto the timescale we're intersted in in
        a=sapply(temp, function(tp){
            tp=as.matrix(tp)
            #print(colnames(propsByTissue))
            #print(colnames(tp))
            #print(tp[1:5,colnames(propsByTissue)])
            #print(j)
            tp[,colnames(propsByTissue)] %*% j
        })
        
    })
})
gene=markerGenes
sim1_v2=lapply(sim1, function(i){
    lapply(i, function(j){
        
        temp=as.data.table(gene)
        temp[, (paste("Sample", 1:length(seq(start.in, stop.in)), sep="_")) := j]

    })
})

save(sim1, file="sim1.RData")
rm(sim1)


####Methods 2 and 3 for simulating RNA-seq datasets:

#The next two methods are based on the mean values from the RNA-seq experiments:
markers <- load.markers(mp) #gets markers that we want. these are currently generated via CIBERSORT
interMarkers=interpolate.markers(markers, boundary.list = bd.list, start = start.in, stop = stop.in) #interpolates onto the timescale we're intersted in in
#numReps=50
noNoise=apply(propsByTissue, 1,  function(j){
  #  lapply(c(1:numReps), function(i){
        lapply(interMarkers, function(tp){
            tp=as.matrix(tp)
            tp[,colnames(propsByTissue)] %*% j
        })
  #  })
})

#################Method 2. Base it on a model of the stdev vs mean relationship
#we can vary the level of noise, but for now, we can do 1
sdlevel=1
meanScale=10
numReps=50
#previously, we developed a model of stdev vs mean in bulk RNA-seq datasets in Arabidopsis (on TPMs)
load("modelStdev.RData")
set.seed(1)
sim2=lapply(noNoise, function(i){ #tissueProps
    lapply(c(1:numReps), function(j){ #numReps
        sapply(i, function(k){ #timepoints
            sapply(k, function(mval){
                mval=mval*meanScale
                sdTemp=sdlevel*max(0.0000001, exp(log(mval)*modelStdev$coefficients[2]+modelStdev$coefficients[1]))
                max(0, mval+rnorm(1, mean=0, sd=sdTemp))
            })
        })
    })
})
save(sim2, file="sim2.RData")
###To make sure this code works, lets calculate the mean/stdev relationship among replicates:
means=sapply(1:length(sim2), function(i){
    sapply(1:length(sim2[[1]][[1]]), function(j){
        mean(sapply(1:length(sim2[[1]]), function(k){
            sim2[[i]][[k]][j]
        }))
    })
})
sds=sapply(1:length(sim2), function(i){
    sapply(1:length(sim2[[1]][[1]]), function(j){
        sd(sapply(1:length(sim2[[1]]), function(k){
            sim2[[i]][[k]][j]
        }))
    })
})
plot(log(means), log(sds), col=rgb(0.5, 0.5, 0.5, 0.1))
save(sim2, file=paste("sim2_part1_sdLevel_", sdlevel, "_meanScale_", meanScale, ".RData", sep=""))




########3. Base it on the covariance between genes in whole seedlings
sdlevel=1
meanScale=10
numReps=50

load("covMat_markers.RData")
set.seed(1)
sim3=lapply(noNoise, function(i){ #tissueProps
   # lapply(i, function(j){ #numReps
        lapply(i, function(k){ #timepoints
            t(mvrnorm(n = numReps, k*meanScale, sdlevel*covMat, empirical = FALSE))
        })
  #  })
})


save(sim3, file=paste("sim3_part1_sdLevel_", sdlevel, "_meanScale_", meanScale, ".RData", sep=""))











