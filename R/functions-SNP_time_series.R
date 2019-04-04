

# checks validity of the input signature for initialize_SNP_time_series method
#
validity_initialize_SNP_time_series=function(chr,pos,base.freq,lib.freqs,
                                             pop.ident, pop.generation,
                                             minfreqchange, minrepl,
                                             max.minor.freq,
                                             winsize,
                                             use.libs,
                                             min.minor.freq,
                                             min.lib.frac,
                                             win.scale,pos.cM)
  
{
  # parameter presence
  if (is.null(chr)) { stop("Parameter 'chr' missing with no default.") }
  if (is.null(pos))  { stop("Parameter 'pos' missing with no default.") }
  if (is.null(base.freq))  { stop("Parameter 'base.freq' missing with no default.") }
  if (is.null(pop.ident))  { stop("Parameter 'pop.ident' missing with no default.") }
  if (is.null(pop.generation)) { stop("Parameter 'pop.generation' missing with no default.") }
  if (is.null(use.libs)) { stop("Parameter 'use.libs' missing with no default.") }
  if (is.null(min.lib.frac)) { stop("Parameter 'min.lib.frac' missing with no default.") }
  if (is.null(win.scale)) { stop("Parameter 'win.scale' missing with no default.") }
  
  if (win.scale=="cM"){
    if (!is.numeric(pos.cM)) stop("'pos.cM' must be a numeric vector.")
    if (length(pos.cM)!=length(chr)) stop("'pos.cM' must have the same length as 'chr', 'pos', 'base.freq' 
                                          and number of rows as in 'lib.freqs'.")
  } else if (win.scale=="bp"){
    pos.cM=rep(NA,length(chr))
  } else stop("'win.scale' must either be 'bp' or 'cM'.")
  
  #----------------
  # test the format of the input variables
  if (!(is.numeric(pos) & is.numeric(base.freq) & is.numeric(pop.ident)& is.numeric(min.minor.freq)
        & is.numeric(max.minor.freq) & is.numeric(minfreqchange) & is.numeric(minrepl)
        & is.numeric(min.lib.frac))) {
    stop("Parameters pos, base.freq, pop.ident, min.minor.freq, max.minor.freq, minfreqchange, 
         min.lib.frac, minrepl have to be numeric.")
  }
  if (!is.logical(use.libs)) stop("Parameter 'use.libs' has to be a logical vector.")
  if (!(min.minor.freq>=0 & min.minor.freq<=1 & max.minor.freq>=0 & max.minor.freq<=1 
        & minfreqchange>=0 & minfreqchange<=1 & min.lib.frac>=0 & min.lib.frac<=1
        & min(base.freq,na.rm = T)>=0 & max(base.freq,na.rm = T)<=1)) { 
    stop("Parameters base.freq, min.minor.freq, max.minor.freq and minfreqchange have to be >=0 and <=1.")
  }
  
  if (!(winsize%%1==0 & minrepl%%1==0)) {
    stop("Parameters winsize, minrepl have to be integers.")
  }
  # the length of the dataset is defined as the number of SNPs and therfore the length of chr, pos
  # , base.freq and the of rows should be equal
  n.SNPs=length(chr)
  if (n.SNPs!=length(pos) | n.SNPs!=length(base.freq) | n.SNPs!=nrow(lib.freqs))
  {
    stop(paste("The length of the vectors chr, pos and base.freq should all be the same as the number of rows in lib.freqs.
               \n length chr:",length(chr),
               "\n length pos:",length(pos),
               "\n length base.freq:",length(base.freq),
               "\n rows in lib.freqs:",nrow(lib.freqs)))
  }
  if (ncol(lib.freqs)!=length(pop.ident) | length(pop.generation)!=length(use.libs) | length(use.libs)!=length(pop.ident))
  {
    stop("Number of libraries in lib.freq, pop.ident and use.libs must all be identical.")
  }
  
  if ((minfreqchange)<0.2 ) {
    warning("The set minfreqchange is below the recommend change
            of 0.2, keep this in mind when working with the results.")
  }
  
  if ((minrepl)<3 ) {
    warning("The set minrepl parameter is below the recommend number
            of 3, keep this in mind when working with the results.")
  }
  
  if ((min.minor.freq)!=0){
    warning("Note that this parameter setting was not included in
            the original publication and should therefore
            only be used for exploratory purposes.")
  }
  
  if ((max.minor.freq-min.minor.freq)>0.3){
    warning("A too large range of starting frequencies is not 
          recommended because it can result in too noisy correlations.
          Consider to modify parameters 'min.minor.freq' and 'max.minor.freq'.
          This parameter setting should only be used for exploratory purposes.")
  } else if ((max.minor.freq-min.minor.freq)>0.05){
    warning("A too large range of starting frequencies is not 
          recommended because it can result in too noisy correlations.
          Consider to modify parameters 'min.minor.freq' and 'max.minor.freq'")
  }
  if ((max.minor.freq<min.minor.freq)){
    stop("Please swap values for 'min.minor.freq' 'max.minor.freq'.")
  }
  
  if (ncol(lib.freqs)<10) {
    warning("Correlations are expected to be unreliable if too few libraries are available.")
  }
  
  pos.cM
}


#' Initialization of time series data as input for haplotype reconstruction
#'
#' This function initializes a genome-wide time series data set that can be used as input for haplotype-block
#' reconstruction.
#'
#' @author Susanne U. Franssen
#' @details The function takes as input genome-wide frequencies of SNPs polarized for the minor frequency allele
#' in the experimnetal starting population for multiple time points and replicates. SNP positions are filtered for a
#' maximum frequency in the experimental starting population and a minimum frequency change in at least one time point
#' for a specified number of replicates. The initialized data is returned as a SNP_time_series object that is
#' required as input for the function \code{reconstruct_hb} to reconstruct unknown haplotype-blocks from the experimental
#' starting population.
#' @export
#' @import foreach data.table
#' @param chr character vector specifying the chromosome name for each genome-wide SNP
#' @param pos numeric vector specifying the chromosomal position for each genome-wide SNP
#' @param base.freq numeric vector specifying the frequency of the minor allele polarized in experimental
#' starting population for each genome-wide SNP
#' @param lib.freqs matrix specifying the frequencies of all genome-wide SNPs (rows) for all different libraries
#' (time points and replicates, columns).
#' @param pos.cM numeric vector corresponding to SNP positions in col.info with
#' genetic positions in cM.
#' @param pop.ident numeric vector specifying the identity of each library in terms of replicate ID
#' @param pop.generation numeric vector specifying the time point of the respective library
#' @param use.libs logical vector specifying which libraries should be used for haplotype-block 
#' reconstruction. The choice taken here determines SNP filtering as parameters \code{minfreqchange}
#' and \code{minrepl} depend on the choice of the data set here. For visualization of marker frequencies,
#' however, the remaining libraries will also be available.
#' @param min.minor.freq numeric specifying the minimum frequency of the minor allele (polarized in the
#' experimental starting population) to be included in the analysis (default=0).
#' @param max.minor.freq numeric specifying the maximum frequency of the minor allele (polarized in the
#' experimental starting population) to be included in the analysis
#' @param winsize numeric specifying the window size on which to perform the analysis
#' @param win.scale character string specifying which genome-wide distance 
#' measure is used for window definition. Options are "bp" (base pairs) or 
#' "cM" (centi Morgan). cM distances can only be used if gentic positions are 
#' provided in 'pos.cM' (default="Mb").
#' @param min.lib.frac minimum fraction of non-NA values for a SNP across 
#' libraries (only using libraries specified in \code{use.libs}) (default=0.75).
#' @param minfreqchange numeric specifying the minimum frequency change required in 'minrepl' replicates required to
#' include the SNP in the analysis
#' @param minrepl numeric specifying the number of replicates, in which the 'minfreqchange' is required to include
#' the SNP in the analysis
#' @return an object of the class \code{SNP_time_series} data
#' @references Franssen, Barton & Schloetterer 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @seealso \code{\link{ex_dat}} \code{\link{SNP_time_series}}
#'
# chr=as.character(example.dat_3$wchr); pos=example.dat_3$gpos; base.freq=example.dat_3$Ball_minfreq;lib.freqs=as.matrix(example.dat_3[,7:55,with=F]);pop.ident=c(1,2,3,4,1,5,2,3,4,1,4,5,2,3,4,1,4,3,4,4,1,5,2,3,4,4,4,4,1,5,2,3,4,4,1,5,2,3,1,5,3,2,4,4,1,4,5,3,2);pop.generation=c(0,0,0,5,7,7,7,7,9,15,15,15,15,15,19,23,25,27,30,35,37,37,37,37,40,45,50,55,59,59,59,59,60,65,75,75,75,75,89,89,89,89,75,89,125,125,125,125,125);use.libs=c(rep(T,44),rep(F,5))
# win.scale="cM";min.minor.freq = 0; max.minor.freq = 3/113; minfreqchange = 0.25; minrepl = 3;pos.cM=1:1682948;min.lib.frac=0.75
initialize_SNP_time_series=function(chr,pos,base.freq,lib.freqs,
                                    pop.ident, pop.generation, use.libs,
                                    minfreqchange=0.2, 
                                    minrepl=3,
                                    max.minor.freq=3/200, 
                                    winsize=500000,
                                    min.minor.freq=0,
                                    min.lib.frac=0.75,
                                    win.scale="bp",pos.cM=NULL)
{
  pos.cM=validity_initialize_SNP_time_series(chr=chr,pos=pos,base.freq=base.freq,lib.freqs=lib.freqs,
                                             pop.ident=pop.ident, pop.generation=pop.generation,
                                             minfreqchange=minfreqchange, minrepl=minrepl,
                                             max.minor.freq=max.minor.freq,
                                             winsize=winsize,
                                             use.libs=use.libs,
                                             min.minor.freq=min.minor.freq,
                                             min.lib.frac=min.lib.frac,
                                             win.scale=win.scale,pos.cM=pos.cM)
#  pos.cM=rep(NA,length(chr))
  
  #-----------------
  # format input variables
  col.info=data.table(chr,pos,base.freq,pos.cM)
  lib.freqs=data.table(lib.freqs)
  nrow(col.info)
  
  # filter for low frequency SNPs in the base (maximum base.freq)
  lib.freqs=lib.freqs[which(col.info$base.freq>=min.minor.freq & col.info$base.freq<=max.minor.freq)]
  col.info=col.info[which(col.info$base.freq>=min.minor.freq & base.freq<=max.minor.freq)]
  nrow(col.info)
  
  # filter for SNPs that increase during the experiment
  # vector counting in how many pops the minfreqchange is exceeded
  pop <- NULL # declared in subsequent foreach 
  pop.counts <- foreach(pop=sort(unique(pop.ident)), .combine="+") %do%
  {
    # returns vector of counts in how many libs for the respective population if minfreqchange is met
    freqchange=foreach(col=lib.freqs[,use.libs & pop.ident==pop,with=FALSE], .combine="+") %do%
    {
      x=(col-col.info$base.freq)>minfreqchange
      x[is.na(x)]=F #change enries with NA to FALSE so that the SNP site is not by default lost
      x
    }
    
    # bolean vector if minfreqchange was reached in at least one library of the respective population
    freqchange=freqchange>0
    # change NAs to FALSE
    #freqchange[is.na(freqchange)]=F # not needed any more as NA handling in previous foreach loop
    if (sum(is.na(freqchange)) > 0) stop("@me check why NAs can still be present in freqchange.")
    freqchange
  }
  lib.freqs=lib.freqs[which(pop.counts>=minrepl)]
  col.info=col.info[which(pop.counts>=minrepl)]
  nrow(col.info)
  
  # filter out SNPs that have to many NAs across libraries
  take=apply(!is.na(lib.freqs),1,sum)/ncol(lib.freqs)> min.lib.frac
  lib.freqs=lib.freqs[take,]
  col.info=col.info[take,]
  nrow(col.info)
  
  # check if a sufficient number of SNPs has remained
  if (nrow(col.info)<1) {
    stop("No SNPs remain after filtering.
         Adjust parameters if possible.")
  } else if (nrow(col.info)<50) {
    warning(paste("After filtering only ",nrow(col.info)," SNPs are left!
                  Maybe check if filtering criteria can be relaxed."))
  } else if (nrow(col.info)<10) {
    warning(paste("After filtering only ",nrow(col.info)," SNPs are left!
                  Rather not continue the analysis with the used parameter setting unless
                  you exactly know what you are doing."))
  }
    
  # sort data
  lib.freqs=lib.freqs[order(col.info$chr,col.info$pos)]
  col.info=col.info[order(col.info$chr,col.info$pos)]
  
  # add window information
  for (c.chr in unique(col.info$chr))
  {
    if (win.scale=="bp") {positions=col.info[chr==c.chr,pos]
    } else {positions=col.info[chr==c.chr,pos.cM]}
    win <- NULL # silly as win is initialized in the subsequent line
    col.info[chr==c.chr,win:=get_win(positions,winsize/2)]
  }
  
  new("SNP_time_series",col.info=col.info,lib.freqs=lib.freqs,pop.ident=pop.ident
      , pop.generation=pop.generation, use.libs=use.libs, winsize=winsize, win.scale=win.scale
      , min.minor.freq=min.minor.freq, max.minor.freq=max.minor.freq
      , minfreqchange=minfreqchange, minrepl=minrepl)
}


# used in initialize_SNP_time_series
#
# returns a vector of window identities for the given 'windowsize' that are non-overlapping and correspond to the provided 'positions'
get_win=function(positions,winsize,offs=0)
{
  breaks=0:ceiling(max(positions)/winsize)*winsize
  win=cut(positions,breaks=breaks+offs,include.lowest=T,labels=F)
  return(win)
}



#------------------------------------------------



# checks validity of the input signature for reconstruct_hb method
#
validity_reconstruct_hb=function(object, chrom, min.cl.size, min.cl.cor, min.inter, single.win)          
{
  # parameter presence
  if (is.null(object)|is.null(chrom)|is.null(min.cl.size)|is.null(min.cl.cor)
      |is.null(min.inter)|is.null(single.win))
  { stop("Parameter missing with no default.") }

  if (!class(object)=="SNP_time_series")
  { stop("object has to be an object of the class SNP_time_series.")
  }
  if (!chrom %in% unique(object@col.info$chr))
  { stop("Value of the parameter chrom must be a chromosome in the SNP_time_series object.")
  }
  if (!(is.numeric(min.cl.size) & is.numeric(min.cl.cor) & is.numeric(min.inter)))
  { stop("Parameters min.cl.size, min.cl.cor and min.inter have to be numeric.")
  }
  
  if (min.cl.size<4) warning("min.cl.size below 4 is not recommended.")
  if (min.inter<2) warning("min.inter below 2 is not recommended.")
  if (min.inter>min.cl.size) stop("min.inter should be <= min.cl.size.")
  if (min.cl.cor<0.4) { 
    #stop("min.cl.cor below 0.4 should not be used! ")
    warning("min.cl.cor below 0.4 should not be used! 
            Be aware that the used correlation threshold is extremely weak
            and the chance is high that unrelateted things are clustered together!")
  } else if (min.cl.cor<0.6) {
    warning("min.cl.cor below 0.6 is not recommended.")
  }
}  




# used in reconstruct_hb
# 
get_clusters=function(dsub.cor,min.cl.size, min.cl.cor, pos.cor)
{
  dsub.cor=dsub.cor[!is.na(dsub.cor[1,]),!is.na(dsub.cor[1,])] # remove NA lines
  
  if (dim(dsub.cor)[1]>1) # only if markers are present that fulfill the criteria
  {
    if (pos.cor) # remove negative corr and replace with 0
    {
      #print("pos.cor")
      dsub.cor[dsub.cor<0]=0 #remove negative correlations
    } #else {print("no pos.cor")}
    
    # getting clusters with a minimum avarage correlation value and a minimum size
    dsub.cordist=1-as.dist(dsub.cor) # convert to dist format
    dsub.clust=hclust(dsub.cordist,method="average") # perform clustering
    #     plot(dsub.clust); abline(h=1-min.cl.cor)
    dsub.cl=cutree(dsub.clust,h = 1-min.cl.cor) # get cluster id based on a max average distance
    #print(dsub.cl)
    dsub.table=table(dsub.cl)
    dsub.clid=as.integer(names(dsub.table))[dsub.table>=min.cl.size] # id of clusters that have more than min.cl.size elements
    dsub.clsnps=list()
    if (length(dsub.clid)>=1) # only if clusters exist that fulfill the critteria
    {
      for (i in 1:length(dsub.clid)) # get list with gnome-wide pos for each cluster
      {
        dsub.clsnps=lappend(dsub.clsnps,as.integer(names(dsub.cl[dsub.cl==dsub.clid[i]])))
      }
      dsub.clsnps=dsub.clsnps[order(as.integer(summary(dsub.clsnps)[,1]),decreasing=T)] # order clusters with largest one first
    }
    return( dsub.clsnps)
    
  } else {
    return(list())
  }
}


# used in reconstruct_hb
# 
get_clusters_dbScan<-function(dsub, min.cl.size, eps) 
{
  xx=dbscan(dsub, eps=eps, minPts = min.cl.size)
  cls=sort(unique(xx$cluster))
  cls=cls[cls!=0] # remove SNP assigned with 0 (which where not included in any cluster)
  
  dsub.clsnps=list() # list with all clusters described by SNP positions
  
  if (length(cls)>0)
  {
    for (i in cls) # get list with genome-wide pos for each cluster
    {
      dsub.clsnps=lappend(dsub.clsnps,as.numeric(row.names(dsub))[xx$cluster==i])
    }
    dsub.clsnps=dsub.clsnps[order(as.integer(summary(dsub.clsnps)[,1]),decreasing=T)] # order clusters with largest one first
    dsub.clsnps
  } else
  {
    list()
  }
  
  
}



# used in reconstruct_hb / get_clusters
# 
# add an object to a list
lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}


# used in reconstruct_hb
#
# get pairs of nodes (clusters) to be combined along the chromosome
# will be merged if two clusters in two overlapping windows share at least min.inter markers
get_edgelist=function(cl.chr,min.inter,single.win) 
{
  el=c() # edgelist not as matrix yet
  for (i.win in 1:(length(cl.chr)-1)) #
  {
    # check if at least one cluster was called in the two overlapping windows
    if (length(cl.chr[[i.win]]) * length(cl.chr[[i.win+1]]) >0)
    {
      for (i in 1:length(cl.chr[[i.win]])) # go through all clusters in the 1st window
      {
        for (j in 1:length(cl.chr[[i.win+1]])) # go through all clusters in the 2nd window
        {
          inter=intersect(cl.chr[[i.win]][[i]],cl.chr[[i.win+1]][[j]])
          if (length(inter)>=min.inter) # if intersects meets the threshold
          {
            #print(paste0(i.win,"_",i,"   ",i.win+1,"_",j))
            x=paste0(i.win,"_",i)
            y=paste0(i.win+1,"_",j)
            el=c(el,x,y)
          }
        }
      }
    }
  }
  if (single.win)
  {
    nodes=get_nodes(cl.chr) # all nodes
    
    # this includes edges from a nodes with no intersections to others to itself
    # --> not needed when the minimum cluster size min.cl.size is sufficiently small (e.g. 4) as then
    # resulting hbrs are not meaningfull anyways
    for (x in setdiff(nodes, unique(el)))
    {
      el=c(el,rep(x,2))
    }
  }

  #print(head(el))
  edge.list=matrix(el,ncol=2,byrow=T)
  return(edge.list)
}



# used in reconstruct_hb 
#
# get markers for each elongated cluster
get_longcluster_markers=function(graph.obj,cl.chr,mode)
{
  # things to recalucalte
  cl.long=get_gwide_cluster_summary(graph.obj)
  graph.vertex=V(graph.obj)$name
  
  cl.long.m=list() # marker SNPs for each long cluster
  for (i in 1:dim(cl.long)[1])
  {
    vertex=cl.long$vertex[i] # starting vertex for each elongated cluster
    el.cl.vertices=graph.vertex[subcomponent(graph.obj,v=vertex,mode="all")] # all vertices of this elongated cluster
    cl.long.m[[i]]=get_win_cluster_marker( el.cl.vertices , cl.chr, mode=mode)
  }
  return(cl.long.m)
}


# used in reconstruct_hb get_longcluster_marker
#
# get summary information of first element and length of genome wide clusters
get_gwide_cluster_summary=function(graph.obj)
{
  graph.cl=clusters(graph.obj)
  graph.vertex=V(graph.obj)$name
  
  vertex=c();len=c()
  for (i in 1:graph.cl$no)
  {
    x=graph.vertex[graph.cl$membership==i]
    #print(x);print(length(x))
    vertex=c(vertex,x[1]) # first vertex in cluster
    len=c(len,length(x)) # length of cluster
  }
  index=1:graph.cl$no
  cl.sum=data.frame(cbind(index,vertex,len))
  cl.sum$vertex=as.character(cl.sum$vertex)
  cl.sum$len=as.integer(as.character(cl.sum$len))
  reorder=order(as.integer(matrix(unlist(strsplit(cl.sum[,2],split = "_")),ncol=2,byrow = T)[,1]))
  cl.sum=cl.sum[reorder,]
  return(cl.sum)
}


# used in reconstruct_hb get_longcluster_marker
#
# returns markers of all given vertices
get_win_cluster_marker=function(vertices,cl.chr,mode="all")
{
  # matrix: 1stcol: window index, 2ndcol: cl index
  x=matrix(as.integer(unlist(strsplit(vertices ,split=("_")))),ncol=2,byrow=T) # col: win, cl
  res=c()
  for (i in 1:dim(x)[1])
  {
    res=c(res,cl.chr[[x[i,1]]][[x[i,2]]])
  }
  
  if (mode=="dupl") # only take markers that were part of the cluster in both windows
  {
    res=res[duplicated(res)]
  } else {
    res=unique(sort(res))
  }
  
  return(res)
}


# used in reconstruct_hb
#
# initialize all nodes of the graph all clusters for each window: X_Y with X the window number and Y the cluster number in the respective window
get_nodes=function(cl.chr)
{
  nodes=c()
  for (i in 1:length(cl.chr))
  {
    if (length(cl.chr[[i]])>0)
    {
      nodes=c(nodes,paste(i,1:length(cl.chr[[i]]),sep="_"))
    }
  }
  return(nodes)
}
