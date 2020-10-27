###generate simulated tissue mixtures:

#start with leaf/root mixtures at 7 days

# grab input data
mix.data <- fread('originalExpressionData/christoff_data.txt')

# scale and normalize the whole thing
for (col.name in setdiff(names(mix.data), "gene")) { #scaling mix testing columns
    set(mix.data, j = col.name, value = scale(mix.data[, col.name, with = F]))
}


leaf <- meta[clean.tissue == 'leaf' & time == '7d']
root <- meta[clean.tissue == 'root' & time == '7d']

percents <- seq(0, 1, 0.1)
perturbations <- seq(0,0.4, .05) 

total.iterations = length(percents) * length(perturbations)
prog <- txtProgressBar(0, total.iterations, style = 3) #create progress bar 
prog.i <- 0
set.seed(1)
count=1
cn=as.character(colnames(mix.data))
leaf.root.test <- rbindlist(lapply(percents, function(mix.x) {
    rbindlist(lapply(perturbations, function(perturb){
        #rbindlist(lapply(1:10, function(count) {
        pairings <- data.table(expand.grid(leaf$sample, root$sample))
        setnames(pairings, c('x', 'y'))
        pairings[, c('tissue_x', 'tissue_y', 'percent_x', 'percent_y') := .('leaf', 'root', mix.x, 1-mix.x)]
        
        pair.out <- pairings[,  runif(length(percent_x), 1 - perturb, 1 + perturb) * percent_x * mix.data[, x, with = F] + runif(length(percent_y), 1 - perturb, 1 + perturb) * percent_y * mix.data[, y, with = F]] #, by = .N]
        
        setnames(pair.out, paste("sample", 1:9, sep = "_"))
        
        pair.out[, gene := mix.data$gene]
        
        #mp <- paste(getwd(), sig.matrix, sep = "/")
        mp <- "signatures/scaled_time_signatures.txt"
        tp <- pair.out
        
        outputs <- main(marker.path = mp, test.samples = tp, bd.list = boundary.list, start.in = 0, return.ratios = TRUE, stop.in = 60, error.fun = 'mono', inf.scale=100)#error.fun = error.inter, inf.scale = scaling)
        outputs[, c('noise', 'mix.x') := .(perturb, mix.x)]
        
        prog.i <<- prog.i + 1
        setTxtProgressBar(prog, prog.i)
        outputs
        # }))  
    }))
}))

prog.i <- 0
close(prog)

