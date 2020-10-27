#######Read in the simulated data and run CoreTissueTimer on it:

simulations=c("sim1.RData")


#lapply(sim1, function)

#example boundary list for some of our data
boundary.list <- list(apex_inflorescence = c(14, 28), apex_vegetative = c(NA,21), flower_whole = c(14, NA), 
                      hypocotyl_hypocotyl = c(NA, 15), leaf_leaf = c(NA,NA), root_root_whole = c(NA,NA), stem_stem = c(NA, NA))



#mp <- "~/Documents/Alan Turing/Projects/data/signatures_from_cibersort.txt"
#tp <- "~/Documents/Alan Turing/Projects/data/christoff_data.txt"
mp <- "signatures/scaled_time_signatures.txt"
#tp <- "originalExpressionData/christoff_data.txt"
#tp <- fread("originalExpressionData/christoff_data.txt")



#tp2 <- as.data.table(sim1[[1]][[1]])
gene=markerGenes
#temp=as.data.table(gene)
#temp[, (paste("Sample", 1:length(seq(start.in, stop.in)), sep="_")) := as.data.table(sim1[[1]][[1]])]

## testing on training data
#outputs <- main(marker.path = mp, test.samples = temp, bd.list = boundary.list, start.in = 0, stop.in = 28, error.fun = 'mono')

outputs_collective_sim1=lapply(sim1, function(i){
    lapply(i, function(j){
        temp=as.data.table(gene)
        temp[, (paste("Sample", 1:length(seq(start.in, stop.in)), sep="_")) := as.data.table(j)]
        main(marker.path = mp, test.samples = temp, bd.list = boundary.list, start.in = 0, stop.in = 14, error.fun = 'mono')
        
    })
})

save(outputs_collective_sim1, file="outputs_collective_sim1.RData")




outputs_collective_sim2=lapply(sim2, function(i){
    lapply(i, function(j){
        temp=as.data.table(gene)
        temp[, (paste("Sample", 1:length(seq(start.in, stop.in)), sep="_")) := as.data.table(j)]
        main(marker.path = mp, test.samples = temp, bd.list = boundary.list, start.in = 0, stop.in = 14, error.fun = 'mono')
        
    })
})

save(outputs_collective_sim2, file="outputs_collective_sim2_rep1to3.RData")


#reorder sim3
sim3_v2=lapply(sim3, function(i){
    lapply(c(1:numReps), function(j){
        sapply(i, function(k){
            k[,j]
        })
    })
})

outputs_collective_sim3=lapply(sim3_v2, function(i){
    lapply(i, function(j){
        temp=as.data.table(gene)
        temp[, (paste("Sample", 1:length(seq(start.in, stop.in)), sep="_")) := as.data.table(j)]
        main(marker.path = mp, test.samples = temp, bd.list = boundary.list, start.in = 0, stop.in = 14, error.fun = 'mono')
        
        
        #temp[, (paste("Sample", 1:length(seq(start.in, stop.in)), sep="_")) := as.data.table(j)]
        #main(marker.path = mp, test.samples = temp, bd.list = boundary.list, start.in = 0, stop.in = 14, error.fun = 'mono')
        
    })
})

save(outputs_collective_sim3, file="outputs_collective_sim3_rep50.RData")




