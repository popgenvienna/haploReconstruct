
####################################
# Generics for SNP_time_series
####################################


setGeneric("reconstruct_hb", function(object, ...) standardGeneric("reconstruct_hb"))


####################################
# Generics for hbr
####################################

setGeneric("plot_cluster_trajectories", function(object, window, ylim=c(0,1), p.median=T, tp.cor=T) {
            standardGeneric("plot_cluster_trajectories") })
setGeneric("plot_marker_trajectories", function(object, hbr_id, ylim=c(0,1), loc.col=T, tp.cor=T) {
  standardGeneric("plot_marker_trajectories") })
setGeneric("map", function(object) standardGeneric("map"))
setGeneric("rev_map", function(object) standardGeneric("rev_map"))
setGeneric("markers", function(object, hbr_id) standardGeneric("markers"))
setGeneric("number_hbr", function(object) standardGeneric("number_hbr"))
setGeneric("plot_hbr_freq", 
           function(object, hbr_id=1, replicate, timepoint, window=1, cols=NULL, add=F
                    , sumstat="mean", cex=0.7
                    , xlab="Genomic position [Mb]", ylab="Marker frequency", xlim=NULL
                    , ylim=c(0,1), pch=20, lwd=2
                    , xaxt=NULL, yaxt=NULL) standardGeneric("plot_hbr_freq"))
setGeneric("inspect_window", function(object, window, colCluster=T) standardGeneric("inspect_window"))
setGeneric("inspect_window_PCA", function(object, window) standardGeneric("inspect_window_PCA"))
setGeneric("inspect_window_avLink", 
           function(object, window, min.cl.cor=0, plotDendro=T, plotCluster=T) standardGeneric("inspect_window_avLink"))
setGeneric("inspect_window_dbScan", 
           function(object, window, minPts=object@min.cl.size, eps, plotkNN=T, plotCluster=T) standardGeneric("inspect_window_dbScan"))


