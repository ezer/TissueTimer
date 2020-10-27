require(data.table)
require(futile.logger)
require(nnls)
require(ggplot2)

## Tissue Timer Main Class
## Aaron Solomon and Daphne Ezer
## Created 5-23-18

#' Imports tissue marker vectors from the file specified in filepath
#' 
#' Expects n x d + 1 matrix in the form of a csv or txt file (filepath), 
#' where n is the number of marker genes
#' and d is the number of tissue types. The extra column must be labeled 
#' "gene" and contain the names of the marker genes.
#' 
#' @param filepath The filepath specifying the location of the csv or txt file for import. 
#' 
#' @return A data.table containing the marker genes and their expression values.
#' 
#' @export

load.markers <- function(filepath = NULL) {
  if (is.null(filepath)) {
    stop("Marker gene filepath was not provided. No genes were imported.")
  }
  
  markers <- fread(filepath) #imports marker genes
  
  if (is.null(markers)) {
    stop("Marker gene file failed to load.")
  } else if (any(dim(markers)) == 0){
    stop("The marker gene file does not provide both tissues and genes.")
  } else if (!("gene" %in% names(markers))){
    stop("No gene column specified.")
  }
  
  return(markers)
}

#' Generates interpolated marker vectors using splines
#' 
#' @param markers A data.table with first column gene and following columns in the format 
#' "tis_sue_nameDATEdate", where date is the number of days of development at which the sample was taken
#' @param boundary.list A named list, where each element of the list is a length 2 vector (c(start,stop))
#' denoting - for the named tissue specified - the first time point at which that tissue begins to grow,
#' and the last time point at which that tissue still appears in any form. the interpolation will
#' ramp up and down from these values. if a value is NA, the tissue will not have the interpolation ramped down to zero, 
#' and instead the most extreme value will be assumed to persist until the edge of the time boundary in that direction
#' if no boundary list is passed, all extreme values are assumed to persist until either bounded edge
#' @param start which day to start on
#' @param stop which day to stop on
#' @param step what interval (in days) to step by 
#' @param scaling.factor a currently aribtrary scaling factor that decides how quickly we drop off tissues at edges 
#' @param ret.recal boolean that determines if we terminate the function early and return recalibrated values - for testing
#' @param plot.recal boolean that determines if we save plots of the recalibrated values - for testing
#' 
#' @return A list of n x d matrices, each of which corresponds to the vector for n genes in d tissues 
#' at the time specified by the list index name.
#' 
#' @export


##
#TO DO: test interpolation. test boundaries, test optimal spline method, 
#DELETE ret.recal param, only there for testing
##
# prototype list for the current tissue: 

interpolate.markers <- function(markers, boundary.list = NULL, start = NULL, stop = NULL, step = 1, scaling.factor = 1, ret.recal = FALSE, plot.recal = TRUE, infinite.scale = 1) {
  
  if (is.null(stop) | is.null(start)) {
    stop("Please provide stop and start times for the sequence.")
  }
  
  tissues <- unique(gsub("DATE[0-9]*$", "", setdiff(unique(names(markers)), "gene"))) #finds unique list of all tissues here
  timepoints <- sort(unique(as.numeric(gsub("[a-zA-Z_]*DATE", "", setdiff(unique(names(markers)), "gene"))))) #sorted list of all timepoints 
  
  inter.points <- seq(start, stop, step) #the points at which we wish to have interpolated values
  
  recalibrated <- rbindlist(lapply(tissues, function(tissue){ #interpolates genes for each tissue and each timepoint
    dt <- markers[, grepl(paste('gene', tissue, sep='|'), names(markers)), with = FALSE]
    times <- as.numeric(gsub("[a-zA-Z_]*DATE", "", setdiff(unique(names(dt)), "gene")))
    setnames(dt, setdiff(unique(names(dt)), "gene"), paste("t", as.character(times),sep=""))
    genes <- markers$gene #removes genes to allow easier matrix computation
    dt[, gene := NULL]
    start.pt <- boundary.list[[tissue]][1]
    stop.pt <- boundary.list[[tissue]][2]
    
    # we're going to change over the version below to incorporate the idea of scaling up the tissue we want to dissapear insteda 
    # of dropping it down -> thus we'll get the actual tissue ratio to go down
    
    scale.factor <- scaling.factor # how much we need to blow up the most extreme tissue to make the scaling work, abitrary val in parameters right now, test this eventually
    
    if (!is.na(start.pt)) { #create 0-bound at left time extreme 
      if (any(start.pt > times)) { #if the user tries to set the 0 column in between available data
        stop(paste("Start point for", tissue, "attempts to set left-extreme expression to zero between available data. Interpolation terminated."))
      } else { # if the user has for some reason specified that the gene expression starts at 0 on a day we actually have data for. if so, we take on step back in time, and begin the interpolation from there
        #update for infiniteTaper: right now, we set the val to scale factor * closest point. what we really want to do
        #is set the furthest point along the interpolation (start or stop) to some very large multiple of the most 
        #extreme point, and then maybe set the point at which this changes to equal the closest point to anchor the 
        #inter polation? 
        if (start.pt %in% times) {
          start.pt <- start.pt - step  
        }
        
        sel.pt <- min(times)
        dt[, paste('t', start.pt,sep = '') := dt[[paste('t', sel.pt, sep = '')]]] #sets our start point within this tissue
        
        print(paste(paste('t', start,sep = ''), infinite.scale, paste('t', sel.pt, sep = ''), min(3*as.numeric(dt[[paste('t', sel.pt, sep = '')]]))))
        dt[, paste('t', start,sep = '') := infinite.scale * as.numeric(dt[[paste('t', sel.pt, sep = '')]])] #infiniteTaper added line!!
        #dt[, paste('t', start,sep = '') :=scale.factor * dt[[paste('t', sel.pt, sep = '')]]] #sets the overall start point for the whole time series to prevent wiggly stuff at the end of the spline
      }
    } else { #if we don't want to taper, simply move left most point to the temporal left start point 
      sel.pt <- min(times)
      dt[, paste('t', min(inter.points), sep = '') := dt[[paste('t', sel.pt, sep = '')]]]
    }
    
    if (!is.na(stop.pt)) { #create 0-bound at left time extreme 
      if (any(stop.pt < times)) { #if the given point is really not the rightmost point
        stop(paste("Stop point for", tissue, "attempts to set right-extreme expression to zero between available data. Interpolation terminated."))
      } else {
        if (stop.pt %in% times) {
          stop.pt <- stop.pt + step  
        }
        
        sel.pt <- max(times)
        dt[, paste('t', stop.pt,sep = '') := dt[[paste('t', sel.pt, sep = '')]]] #sets our start point within this tissue
        print(paste("place 2", infinite.scale, paste('t', sel.pt, sep = ''), min(dt[[paste('t', sel.pt, sep = '')]])))
        
        dt[, paste('t', stop,sep = '') := infinite.scale * as.numeric(dt[[paste('t', sel.pt, sep = '')]])]
        #dt[, paste('t', stop,sep = '') := scale.factor * dt[[paste('t', sel.pt, sep = '')]]] #sets the overall start point for the whole time series to prevent wiggly stuff at the end of the spline
      }
    } else { #if we don't want to taper, simply move left most point to the temporal left start point 
      sel.pt <- max(times)
      dt[, paste('t', max(inter.points), sep = '') := dt[[paste('t', sel.pt, sep = '')]]]
    }
    
    inter.data <- rbindlist(lapply(seq(1,nrow(dt), 1), function(row){ #interpolates each gene within the current tissue
      sp <- interpolate(x = gsub("t", "", names(dt)), y = dt[row], tp = inter.points, interpolator = "linear", method = "linear")
      ret <- as.list(sp$y)
      names(ret) <- sp$x
      return(ret)
    }))
    
    inter.data[inter.data < 0] <- 0 
    inter.data[, c('gene','tissue') := .(genes, tissue)]
    
    
    inter.data
  }))
  
  melt.recal <- melt(recalibrated, id.vars = c('gene', 'tissue')) #melts all of the interpolated data into one frame
  #maybe we can just fill bounds with zeros here? 
  setnames(melt.recal, c("variable", "value"), c("time", "tpm"))
  
  plot.recal=F
  if (plot.recal) { # testing method, delete for production
    
    ## output trajectory of interpolated genes to see how we're constructing and fucking up our trajs
    # 
  #   # clusters <- rbindlist(lapply(unique(recalibrated$tissue), function(tis) {
  #   #   
  #   #   curr <- recalibrated[tissue == tis]
  #   #   mat <- as.matrix(curr[, tissue := NULL])
  #   #   rownames(mat) <- mat[, colnames(mat) == "gene"]
  #   #   mat <- mat[, colnames(mat) != "gene"]
  #   #   
  #   #  # mat <-as.numeric(mat)
  #   #  # print(dim(mat))
  #   #   mat=apply(mat, c(1,2), function(temp){as.numeric(temp)})
  #   #   print(dim(mat))
  #   #   print(mat[1:3, 1:3])
  #   #   km <- kmeans(mat, centers = 3)
  #   #   clus <- km$cluster
  #   #   
  #   #   data.table(tissue = tis, gene = names(clus), cluster = clus)
  #   # }))
  #   # 
  #   #clustered.genes <- melt.recal[clusters, on = c("tissue", "gene")]
  #   
  #   # ggplot(data = clustered.genes, aes(x = jitter(as.numeric(time)), y = tpm, color = factor(cluster))) + geom_smooth() + facet_wrap(~tissue) + ggtitle("K-Means Clustered Interpolated Marker Genes") + xlab("Time (days)") + ylab("TPM")
  #   # ggsave(filename = 'curves.png', device = 'png')
  #   # ggplot(data = clustered.genes[time %in% c(0,7,14,21,28)], aes(x = tissue, y = tpm, color = factor(cluster)))+ geom_boxplot() + facet_wrap(~time) + theme(axis.text.x = element_text(angle = 90))
  #   # ggsave(filename = 'boxes.png', device = 'png')
  #   # 
  #   # raw plots
  #   
  #   dt <- melt(markers, id.vars = 'gene')
  #   dt[, c('tissue', 'time') := tstrsplit(dt$variable, split = "DATE")]
  #   dt[, variable := NULL]
  #   dt[, time := factor(time, levels = seq(1,max(as.numeric(time)),1))]
  #   setnames(dt, "value", "tpm")
  #   
  #  # dt <- dt[clusters, on = c("tissue", "gene")]
  #   
  #   #
  #   # ggplot(data = dt, aes(x = time, y = tpm, color = factor(cluster))) +
  #   #   geom_boxplot() + facet_wrap(~tissue) + 
  #   #   ggtitle("K-Means Clustered Raw Marker Genes") + 
  #   #   xlab("Time (days)") + ylab("TPM")
  #   # ggsave(filename = 'raw_boxes.png', device = 'png')
  #   # 
  #   
  #   ## boxes and curves
  #   
  #   #clustered.genes[, time := factor(time, levels = seq(0,max(as.numeric(time)),1))]
  #   
  #   #curves and boxes
  #  # curves_and_boxes <<- ggplot(data = clustered.genes, 
  # #         aes(x = as.numeric(time), y = tpm, color = factor(cluster))) + 
  #  #   geom_smooth() + 
  #   #  geom_boxplot(data = dt, aes(x = time, y = tpm, color = factor(cluster))) + 
  #   #  facet_wrap(~tissue) + 
  #   #  ggtitle("K-Means Clustered Marker Genes") + 
  #   #  xlab("Time (days)") + ylab("TPM") + scale_x_discrete(limits = seq(1,35,1))
  # #  curves_and_boxes
  #   #ggsave(filename = 'boxes_and_curves.png', device = 'png')
  #   
  #   #points
  #  # point_hist <<- ggplot(data = clustered.genes, 
  # #         aes(x = as.numeric(time), y = tpm, color = factor(cluster))) + geom_point(aes(alpha = .2)) + facet_wrap(~tissue) + 
  # #    ggtitle("K-Means Clustered Marker Genes") + 
  # #    xlab("Time (days)") + ylab("TPM") 
  # #  point_hist
  #   #ggsave(filename = 'hist.png', device = 'png')
  #   ## end traj output
  # }
  # 
  }
  export.timings <- sort(unique(melt.recal$time))
  
  #recasts interpolated data so that each frame contains all tissues at a certain time point
  timed.frames <- lapply(export.timings, function(curr.time){
    time.frame <- dcast(data = melt.recal, formula = gene ~ tissue, subset = .(time == curr.time), value.var = "tpm")
    setorder(time.frame, "gene")
    time.frame <- as.matrix(time.frame)
    rownames(time.frame) <- time.frame[,1]
    time.frame <- time.frame[,2:ncol(time.frame)]
    
    genes <- rownames(time.frame)
    time.frame <- apply(time.frame,2, as.numeric)
    rownames(time.frame) <- genes
    
    time.frame
  })
  
  #labels timed frames for return
  names(timed.frames) <- export.timings
  
  return(timed.frames)
}


#' Interpolates points x and y into a function
#' 
#' @param x A series of points x that denote the times over which we have data
#' @param y A series of points y that denote the values at each time 
#' @param tp The time points at which we wish to 
#' @param interpolator fill this out 
#' @param method fill this out 
#' 
#' @return A function that returns interpolated values from the points and times provided. 
#' 
#' @export
#Interpolation function
#returns list with numeric vectors x and y 
interpolate <- function(x, y, tp, interpolator = c('linear', 'spline'), method = c('linear', 'natural')){
  
  if (interpolator == "linear") {
    return(approx(x = x, y = y, xout = tp, method = method, rule = 2))
  } else if (interpolator == "spline") {
    return(spline(x = x, y = y, xout = tp, method = method))
  } else {
    flog.error("No splining function found.")
    stop()
  }
  
}


#todo add pseudo count somewhere

#' Takes sample data and computes most likely tissue ratio
#' 
#' @param sample sample data: a matrix with rownames == genes and one column that is the data
#' @param marker marker data:  n x d matrix of n genes and d tissues.
#' @param time flesh this out 
#' @param name flesh this out 
#' 
#' @return list with ratios (from regression coiefficients) and L2 normed residuals returned
#' 
#' @export
fit.sample <- function(sample, marker, time = NULL, name = NULL) {
  
  prep.start <- proc.time()
  
  marker <- marker[, apply(marker, 2, sd) > 0.01] #cannot regress on constant columns
  
  #filter sample and marker so they have the same genes
  marker <- marker[rownames(marker) %in% rownames(sample),] #need same genes 
  sample <- sample[rownames(marker),]
  
  prep.stop <- proc.time()
  
  prep.time <<- prep.time + prep.stop['elapsed'] - prep.start['elapsed']
  
  nnls.start <- proc.time()
  
  fit <- nnls(marker, sample)
  
  nnls.stop <- proc.time()
  
  nnls.time <<- nnls.time + nnls.stop['elapsed'] - nnls.start['elapsed']
  
  clean.start <- proc.time()
  
  coef <- coefficients(fit) #calciulating and normalizing tissue ratio
  ratio <- coef / sum(coef)
  names(ratio) <- colnames(marker)
  
  resid <- fit$residuals^2 #calculating residuals with 2-norm
  
  clean.stop <- proc.time()
  clean.time <<- clean.time + clean.stop['elapsed'] - clean.start['elapsed']
  
  
  list(ratios = ratio, residuals = resid)
}

#' Runs Tissue Timer. 
#' 
#' @param marker.path a filepath to a .txt or .csv file with format of n x d + 1, 
#' where n is the number of marker genes (rows)
#' and d is the number of tissue types (columns). The extra column must be labeled 
#' "gene" and contain the names of the marker genes. 
#' @param test.samples samples to test the markers on -> specify exactly the format
#' @param bd.list named list for each tissue that specifies whether
#' the tissue dissapears from the planet and should be tapered or whether it is 
#' persistent and its expression values should be allowed to remain indefinitely 
#' @param start.in when to allow the interpolation to begin
#' @param stop.in when to stop interpolation
#' @param sc.factor arbitrary tapering factor, will test and adjust 
# 
# 
#' @export

main <- function(marker.path = NULL, test.samples = NULL,
                 bd.list = list(apex_inflorescence = c(14, 28), apex_vegetative = c(NA,21), 
                                flower_whole = c(21, NA), hypocotyl_hypocotyl = c(NA, 15), 
                                leaf_leaf = c(NA,NA), root_root_whole = c(NA,NA), 
                                stem_stem = c(NA, NA)), 
                 start.in = NULL, stop.in = NULL, step.in = 1, sc.factor = 4, return.ratios = FALSE, error.fun = 'spline', inf.scale = 1000) {
  
  fit.time <- 0
  spline.time <- 0
  op.time <- 0
  prep.time <<- 0
  nnls.time <<- 0
  clean.time <<- 0
  
  
  
  # flog.info("Main process initialization")
  initialize.time <- proc.time()
  
  markers <- load.markers(marker.path) #gets markers that we want. these are currently generated via CIBERSORT
  
  # flog.info("Markers loaded succesfully")
  
  if (is.data.table(test.samples)) {
    test.data <<- test.samples
  } else {
    test.data <<- fread(test.samples)
  }
  
  for (col.name in setdiff(names(test.data), "gene")) { #scaling test data columns
    set(test.data, j = col.name, value = scale(test.data[, col.name, with = F]))
  }
  
  # flog.info("Test data loaded succesfully")
  
  inter.marks <<- interpolate.markers(markers, boundary.list = bd.list, start = start.in, stop = stop.in, step = step.in, scaling.factor = sc.factor, infinite.scale = inf.scale) #interpolates onto the timescale we're intersted in in
  #saveRDS(object = inter.marks, file = "~/Documents/Alan Turing/Projects/TissueTimerView/data/intermarks.RDS")
  
  inter.time <- proc.time()
  # flog.info(paste("Interpolation complete - time elapsed:", round(inter.time['elapsed'] - initialize.time['elapsed'], digits = 3)))
  
  outputs <- rbindlist(lapply(setdiff(names(test.data), 'gene'), function(samp.name){ #for each sample in the testing data
    
    #for each time point, computes coefficients of fits and residuals, returns
    start.fit <- proc.time()
    fits <- lapply(names(inter.marks), function(timepoint) {
      sample <- as.matrix(test.data[, samp.name, with = F])
      rownames(sample) <- test.data$gene
      fit.sample(sample, marker = inter.marks[[timepoint]], 
                 time = inter.marks[[timepoint]], 
                 name = timepoint)
      
    })
    stop.fit <- proc.time()
    fit.time <<- fit.time + stop.fit['elapsed'] - start.fit['elapsed']
    
    
    names(fits) <- names(inter.marks) #adds names back in after lapply
    resid <- data.table(time = as.numeric(names(fits)), residual = unlist(lapply(fits, function(x){sum(x$residuals)}))) #extracts residuals into data.table -> L2 norm already computed in fit functon so we just sum here
    
    spline.start <- proc.time()
    #
    
    
    if (error.fun == "linear") {
      sp.fun <- approxfun(resid$time, resid$residual, method = "linear") #testing! does linear residualing remove blip points
    } else if (error.fun == "mono") {
      sp.fun <- splinefun(resid$time, resid$residual, method = "monoH.FC") #calculates spline function for use in optimize
    } else {
      sp.fun <- splinefun(resid$time, resid$residual, method = "natural") #calculates spline function for use in optimize
    }
    
    
    spline.stop <- proc.time()
    
    spline.time <<- spline.time + spline.stop['elapsed'] - spline.start['elapsed']
    
    tdf <- data.table(times = seq(start.in, stop.in, .5), vals = sp.fun(seq(start.in, stop.in, .5)))
    
    plot.spline(samp.name, test.data, start.in, stop.in, error.fun)
    
    #ggplot(data=tdf, aes(x = times, y = vals)) + geom_point(aes(alpha = 0.3)) + geom_smooth() + geom_vline(aes(alpha = 0.3), xintercept = as.numeric(keep.d[sample == samp.name, time]))
    #ggsave(filename = paste(samp.name, ".png", sep = "_"), device = "png")
    
    #not officially optimizing here, just using min value from tdf (interpolated)
    # op.start <- proc.time()
    # op.point <- optimize(sp.fun, c(start.in, stop.in)) #calculates minimum point - specify how far we allow it togo
    # op.stop <- proc.time()
    # 
    # op.time <<- op.time + op.stop['elapsed'] - op.start['elapsed']
    # 
    # start.y <- sp.fun(start.in) #also check endpoints for function
    # stop.y <- sp.fun(stop.in)
    # 
    #min.y.idx <- which.min(c(start.y, op.point$objective, stop.y))
    #min.x <- c(start.in, op.point$minimum, stop.in)[min.y.idx]
    
    #below are the min values for linear interpolation.
    #modified again below them to account for multiple points being equal
    #min.x <- tdf[which.min(tdf[, vals]), times]
    #min.y.val <- tdf[which.min(tdf[, vals]), vals]
    
    #here's the one that picks middle points
    min.y <- min(tdf[,vals])
    min.x <- round(mean(tdf[min.y == vals, times])) #this predisposes to be lower by 
    min.y.val <- tdf[times == min.x, vals]
    
    #using the minimized residual time point, return tissue ratio
    ret.point <- names(fits)[which.min(abs(as.numeric(names(fits)) - min.x))] #finds closest timepoint with calculated ratios
    
    rat <- unlist(fits[[ret.point]]['ratios'])
    names(rat) <- gsub("ratios.", "", names(rat), fixed = T)
    
    list(sample.name = rep(samp.name, length.out = length(rat)), 
         time = rep(ret.point, length.out = length(rat)), 
         tissue = names(rat), 
         percent = rat)
  }))
  
  # flog.info(paste("Algorithm complete - time elapsed:", round(proc.time()['elapsed'] - initialize.time['elapsed'], digits = 3)))
  # 
  # flog.info("Runtime Breakdown:")
  # flog.info(paste("\tFit Time:", fit.time))
  # flog.info(paste("\t\tPrep Fit Time:", prep.time))
  # flog.info(paste("\t\tNNLS Time:", nnls.time))
  # flog.info(paste("\t\tClean Fit Time:", clean.time))
  # flog.info(paste("\tSplining Time:", spline.time))
  # flog.info(paste("\tOptimization Time:", op.time))

  
  
  return(outputs)
}

plot.spline <- function(samp.name, test.data, start.in, stop.in, error.fun) {
  
  #inter.marks <- readRDS("~/Documents/Alan Turing/Projects/TissueTimerView/data/intermarks.RDS")
  
  fits <- lapply(names(inter.marks), function(timepoint) {
    sample <- as.matrix(test.data[, samp.name, with = F])
    rownames(sample) <- test.data$gene
    fit.sample(sample, marker = inter.marks[[timepoint]], 
               time = inter.marks[[timepoint]], 
               name = timepoint)
    
  })

  names(fits) <- names(inter.marks) #adds names back in after lapply
  resid <- data.table(time = as.numeric(names(fits)), residual = unlist(lapply(fits, function(x){sum(x$residuals)}))) #extracts residuals into data.table -> L2 norm already computed in fit functon so we just sum here

  if (error.fun == "linear") {
    sp.fun <- approxfun(resid$time, resid$residual, method = "linear") #testing! does linear residualing remove blip points
  } else if (error.fun == "mono") {
    sp.fun <- splinefun(resid$time, resid$residual, method = "monoH.FC") #calculates spline function for use in optimize
  } else {
    sp.fun <- splinefun(resid$time, resid$residual, method = "natural") #calculates spline function for use in optimize
  }
  
  real.error <- data.table(times = resid$time, vals = resid$residual)
  tdf <- data.table(times = seq(start.in, stop.in, .5), vals = sp.fun(seq(start.in, stop.in, .5)))
  
  real.error[, c('type','size', 'alpha') := .('Point', 4, .9)]
  tdf[,  c('type','size', 'alpha') := .('Interpolated', 2.5, 1)]
  
  tdf <- rbind(real.error, tdf)
  
 # ggplot(data=tdf, aes(x = times, y = vals, shape = factor(type), size = size, alpha = alpha)) + 
 #    geom_point(aes()) + geom_vline(aes(alpha = 0.3), xintercept = as.numeric(keep.d[sample == samp.name, time]))
 #  ggsave(filename = paste0(samp.name, "_hyper_tuned.png"), device = "png")
}
