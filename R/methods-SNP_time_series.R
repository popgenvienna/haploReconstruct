

#' Display of a SNP_time_series object
#'
#' This function displays the summarized features of a SNP_time_series object
#'
#' @aliases show.SNP_time_series
#' @author Susanne U. Franssen
#' @details This function displays the summarized features of a SNP_time_series object
#' @export
#' @param object object of the class SNP_time_series data
#' @references Franssen et al. 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @seealso \code{\link{ex_dat}} \code{\link{initialize_SNP_time_series}} \code{\link{SNP_time_series}}
#' 
setMethod("show","SNP_time_series",
            function(object) {
              cat("Object of class",class(object),"\n")
              cat("Number of populations included:",length(unique(object@pop.ident)),"\n")
              cat("Number of libraries included:",length(object@pop.ident),"\n")
              cat("Number of libraries used for filtering:",sum(object@use.libs),"\n")
              cat("Data has been filtered for SNPs with:\n"
                  ,"\t allele frequencies >=",object@min.minor.freq," in the founder population\n"
                  ,"\t allele frequencies <=",object@max.minor.freq," in the founder population\n"
                  ,"\t a frequency increase of",object@minfreqchange," in at least ",object@minrepl,"replicates\n")
              cat("Number of SNPs after filtering:",nrow(object@col.info),"\n")
              cat("Window definition:\n")
              cat("\t distance measure used:",object@win.scale,"\n")
              cat("\t applied windowsize:",object@winsize,"\n")
            })





##object=obj.SNP.ts; chrom="2L";  min.cl.size=4; min.cl.cor=0.6; min.inter=4
#
#' Reconstruction of haplotype-blocks
#'
#' This function reconstructs haplotype-blocks from a SNP_time_series object based on correlated allele
#' frequency changes aacross time points and replicates
#'
#' @aliases reconstruct_hb reconstruct_hb.SNP_time_series
#' @author Susanne U. Franssen
#' @details This function reconstructs haplotype-blocks via markers that show correlated frequency changes across
#' multiple time-points and replicates. In a sliding window based approach SNPs are clustered by average linkage
#' clustering based on correlations of their frequencies across time points and replicates. Clusters are build for
#' each window using the parameters \code{min.cl.size} and \code{min.cl.cor}. Clustered in overlaping windows (window size/2) are
#' elongated to haplotype-blocks if two adjacent clusters have \code{min.inter} identical markers./n
#' Reconstructed haplotype blocks were constructed from at least two overlapping windows and consist of markers that
#' are present in two overlaping window clusters.
#' @export
#' @importFrom stats hclust cor cutree
#' @importFrom methods new
#' @import igraph
#' @param object object of the class SNP_time_series data
#' @param chrom a string specifying the name of the chromosome present in the object for which
#'  haplotype-block should be reconstructed
#' @param min.cl.size numeric specifying the minimum number of correlated markers in a window for 
#' haplotype-block reconstruction
#' @param min.cl.cor numeric specifying the correlation between marker SNPs required using average
#' linkage clustering for markers to be assembled in one cluster
#' @param min.inter numeric specifying the minimum number of markers in two clusters of overlapping 
#' windows required to be identical for cluster elongation across windows to build haplotype-blocks
#' @param single.win boolean specifying that are supported by only a cluster in one window are 
#' also included in the markers of reconstructed blocks. If FALSE only hbrs are returned that span
#' at least two windows and only markers being present in the intersection between overlapping
#' windows are included. 
#' @param transf boolean indicating if frequency data should be square root transformed 
#' prior to calaculation pariwise correlations with Pearson's correlation coefficient.
#' (if transformations are chosen they will be done in the order: sqrt, arcsine, scale)
#' @param arcsine boolean indicating if frequency data should be arcsine transformed 
#' prior to calaculation pariwise correlations with Pearson's correlation coefficient.
#' (if transformations are chosen they will be done in the order: sqrt, arcsine, scale)
#' @param scaleSNP boolean indicating whether time series allele frequency data
#' for each SNP should be scaled (scaling to a mean of zero and standard deviation of 1.
#' (if transformations are chosen they will be done in the order: sqrt, arcsine, scale)
#' @param pos.cor boolean indicating if negative correlations should be set to zero
#' prior to clustering.
#' @param clusterM indicating the clustering method, choose from: "avLink" (average linkage clustering)
#' and "\href{https://CRAN.R-project.org/package=dbscan}{dbscan}".
#' @param eps size of the epsilon neighborhood when clustering with \href{https://CRAN.R-project.org/package=dbscan}{dbscan}.
#' For details please refer to the dbscan package. 
#' @return an object of the class \code{\link{hbr}}
#' @references Franssen et al. 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @seealso \code{\link{ex_dat}} \code{\link{initialize_SNP_time_series}} \code{\link{SNP_time_series}}
#'
setMethod("reconstruct_hb","SNP_time_series",
          function(object, chrom, min.cl.size=4, min.cl.cor=0.8, min.inter=2, single.win=F
                   , transf=T, arcsine=F, scaleSNP=T, pos.cor=F, clusterM="avLink", eps=0.3) {
  
  validity_reconstruct_hb(object, chrom, min.cl.size, min.cl.cor, min.inter, single.win)          
  
  chr <- NULL # sily as in object@col.info[chr==chrom] chr is already defined
  #------------------------
  # calculate clusters for windows
  #
  cl.chr=list() # every list element can be:
  #1) a filled list (each element within represents a cluster always as a vector of genomic 
  #   SNP positions within)
  #2) an empty list: list()
  # a) if SNPs passed the initial filtering criteria (sufficient rise in multiple replicates) 
  #   but cluster criteria (enough markers per cluster) were not fulfilled
  #
  # start and end windows for the respective chromosome
  # win is half the chosen window size and 
  #   therefore two at a time will be considered
  win <- NULL # silly as win is defined in object@col.info
  swin=1; ewin=max(object@col.info[chr==chrom]$win)-1 
  #
  for (i.win in swin:ewin) # go through each window (each made out of 2 subwindows)
  {
    #print(paste("Processing window:",i.win))
    a=NULL
    
    # time series data for a window of "winsize" (equals two windows encoded as win)
    dsub.pos=object@col.info[chr==chrom & (win==i.win | win==i.win+1)]$pos
    dsub=object@lib.freqs[object@col.info$chr==chrom & (object@col.info$win==i.win 
                            | object@col.info$win==i.win+1),object@use.libs,with=F]
    
    if (dim(dsub)[1]>=min.cl.size)
    {
      # calculation of the correlation matrix (correlations between SNPs)
      a=t(dsub)
      colnames(a)=dsub.pos #set col names to genomic positions
      
      if (transf) {a=sqrt(a)}#; print("trans")}
      #print(a[1,1])
      if (arcsine) {a=asin(a)}#; print("asin")}
      #print(a[1,1])
      if (scaleSNP) {a=scale(a)}#; print("scale")}
      #print("--")
      
      if(clusterM=="avLink")
      {
        #print("avLink")
        dsub.cor=cor(a,use="pairwise.complete.obs")
        
        # identifying clusters from the correlation matrix
        dsub.clsnps=get_clusters(dsub.cor,min.cl.size=min.cl.size, min.cl.cor=min.cl.cor, pos.cor=pos.cor)
      } else if (clusterM=="dbScan")
      {
        #print("dbscan")
        dsub.clsnps=get_clusters_dbScan(t(a), min.cl.size=min.cl.size, eps=eps)
      } else {
        stop("The clusterM parameter has to be either avLink or dbScan.")
      }
      
      
      # add clusters for this window
      cl.chr[[length(cl.chr)+1]] <- dsub.clsnps # save all clusters for chromosome-wide windows
    } else {
      cl.chr[[length(cl.chr)+1]] <- list() # add empty list for that window so windows with no identified clusters can be retrieved
    }
  }
  #cl.chr
  
  if (sum(unlist(lapply(cl.chr,length))) > 0) # at least one window contains a cluster
  {
    #------------------------
    # elongate clusters across windows to generate reconstructed haplotype-blocks
    # for cluster prolongation / building graph
    #
    # rows indicate edges between two neighbouring vertices
    edge.list=get_edgelist(cl.chr,min.inter=min.inter,single.win=single.win) 
    graph.obj=graph.edgelist(edge.list, directed=TRUE)
    if (single.win){mode="all"}else{mode="dupl"}
    cl.long.m=get_longcluster_markers(graph.obj,cl.chr,mode=mode) # contains markers of elongated clusters
    
    SNPts=new("SNP_time_series", col.info=object@col.info[chr==chrom]
              , lib.freqs=object@lib.freqs[object@col.info$chr==chrom]
              , pop.ident=object@pop.ident
              , pop.generation=object@pop.generation
              , use.libs=object@use.libs
              , winsize=object@winsize
              , min.minor.freq=object@min.minor.freq
              , max.minor.freq=object@max.minor.freq
              , minfreqchange=object@minfreqchange
              , minrepl=object@minrepl)
    
    # return
    new("hbr",  dat=SNPts
        , chromosome=chrom
        , min.cl.size=min.cl.size, min.cl.cor=min.cl.cor, min.inter=min.inter
        , single.win=single.win
        , transf=transf, arcsine=arcsine, scaleSNP=scaleSNP, pos.cor=pos.cor
        , cl.chr=cl.chr, cl.long.m=cl.long.m)
  } else # no clusters at all were reconstructed
  {
    SNPts=new("SNP_time_series", col.info=object@col.info[chr==chrom]
              , lib.freqs=object@lib.freqs[object@col.info$chr==chrom]
              , pop.ident=object@pop.ident
              , pop.generation=object@pop.generation
              , use.libs=object@use.libs
              , winsize=object@winsize
              , min.minor.freq=object@min.minor.freq
              , max.minor.freq=object@max.minor.freq
              , minfreqchange=object@minfreqchange
              , minrepl=object@minrepl)
    

    # return
    new("hbr",  dat=SNPts
        , chromosome=chrom
        , min.cl.size=min.cl.size, min.cl.cor=min.cl.cor, min.inter=min.inter
        , single.win=single.win
        , transf=transf, arcsine=arcsine, scaleSNP=scaleSNP, pos.cor=pos.cor
        , cl.chr=cl.chr, cl.long.m=list())
  }
  

})



