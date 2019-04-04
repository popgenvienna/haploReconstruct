

#' Data input from a sync file 
#'
#' Reads in SNP time series data from a file with \code{.sync} format. 
#'
#' @author Susanne U. Franssen
#' @details Time series data from a file with \code{sync} format are read in. The sync
#' format is specified in Kofler et al. 2011 (PoPoolation2: identifying differentiation 
#' between populations using sequencing of pooled DNA samples (Pool-Seq)). Allele counts 
#' are read in for each library and SNP and transformed to allele frequencies. Allele 
#' frequencies are polarized for the minor and major allele of a specifies (sub-)set of 
#' libraries, i.e. libraries of the experimentla founder population. Frequencies are 
#' determined only based on the counts of the two most common alleles in the specified 
#' base populations \code{base.pops}.
#' Please note: This procedure does not substitute a proper SNP calling. Provided sync 
#' files are expected only to contain positions of previously called SNPs and at least 
#' two alleles should be present in the specified base populations.
#' @export
#' @import foreach data.table
#' @importFrom stringi stri_split_fixed
#' @importFrom matrixStats rowRanks
#' @param file the name of the ".sync" file where the data should be read from. Sync 
#' files are specified in Kofler et al. (2011). Sync files contain 3 + n columns with; 
#' col 1: chromosome (reference contig), col 2: position (in the reference contig), 
#' col 3: reference allele, col >3: sync entries for allele frequencies for all populations 
#' in the form A-count:T-count:C-count:G-count:N-count:deletion-count.
#' Sync files originally don't have a header but headers are accepted when specified 
#' with \code{header=T}.
#' @param base.pops logical vector with the same length as the number of libraries 
#' present in the sync file. Libraries indicated with TRUE will be used for identification 
#' on the two main alleles (minor and major allele). Allele frequencies of all libraries 
#' will subsequently be polarized for the minor allele in this specified subset.
#' @param mincov minimum coverage to calculate allele frequencies. If the sum of allele 
#' counts of the minor and major allele are below this threshold the respective frequency 
#' will be encoded as NA (default=15).
#' @param header logical value specifying whether a header is present in the provided 
#' sync file.
#' @param polaRise numberical vector specifying polarisation of SNPs based on the rising allele
#' based on the mean frequency change between the given column pairs, e.g. list(c(1,4),c(5,8),c(9,12))
#' the mean frequency change is calclulated between the mean frequency change from cloumn 1 to 4
#' 5 to 8 and 9 to 12 (column count starting from the first library provided). 
#' Polarisation is performed so that the mean frequency change is positive.
#' default: ignore and polarise on the minor allele in the base population
#' @return a data.table with 6 plus N columns with; col 1: chr (chromosome), col 2: pos 
#' (position on respective chromosome), col 3: ref (reference allele), col 4: minallele 
#' (minor allele across all specified base populations), col 5: majallele (major allele 
#' across all specified base populations), col 6: weighted mean frequency of all specified 
#' base populations poloarlized for the minor allele, col >6: allele frequency of the minor 
#' allele for each library
#' @references Franssen, Barton & Schloetterer 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#'
# file="/Volumes/cluster/Suse/haplotype_simul_number/MS03_haplo_reconstruct/01_simulate/test.sync"; base.pops=c(T,T,rep(F,5)); header=T; mincov=15
sync_to_frequencies<-function(file, base.pops, header, mincov=15, polaRise=NULL)
{
  cat("Reading in sync file.")
  x=validity_sync_to_frequencies(file, base.pops, header, mincov)
  # set colnames
  colnames(x)[1:3]<-c("chr","pos","ref")
  if (!header) colnames(x)[4:ncol(x)]<-paste0("L",1:(ncol(x)-3))
  
  # only for my case where alle ref alleles were encoded as T and therefore read in as TRUE
  if (is.logical(x$ref)) x$ref="T" 
  
  # for testing
  #x=x[1:5,1:10]
  
  # list with 4 col (ATCG) matrix for each lib
  x.counts.4alleles=getATCG(x, cols=4:ncol(x))
  
  # create sync entry of all base pops to be fused to get min and maj alleles
  dd <- NULL # dd initialized in subsequent foreach
  base.sum.4alleles <- foreach(dd=x.counts.4alleles[base.pops], .combine="+") %do% dd
  
  # index (1:4 in sync entry) for the minor and major allele for each SNP
  x.allele.index=top2(base.sum.4alleles)
  
  # add summed up counts of all base.pops to the overl list with 4 allele counts for all pops
  all.counts.4alleles=append(list(base.sum.4alleles), x.counts.4alleles)
  
  # list with 2 col (min max allele) for each lib
  d <- NULL # initialized in subsequent foreach
  all.counts.2alleles <- foreach(d=all.counts.4alleles) %do% {
    cbind(d[cbind(1:nrow(x.allele.index), x.allele.index$Min)], d[cbind(1:nrow(x.allele.index), x.allele.index$Max)])
  }
  #names(all.counts.2alleles) <- c("basePops",colnames(x)[4:ncol(x)]) # add lib names
  
  # get nt info for minor and major allele
  minallele <- c("A", "T", "C", "G")[x.allele.index$Min]
  majallele <- c("A", "T", "C", "G")[x.allele.index$Max]
  
  # minfreq for each replicate and timepoint
  lib <- NULL # initialized in subsequent foreach
  minfreqs=matrix(unlist(foreach(lib=all.counts.2alleles) %do% (lib[,1]/(lib[,1]+lib[,2]))), ncol=length(all.counts.2alleles)) 
  colnames(minfreqs)<-c("basePops",colnames(x)[4:ncol(x)]) # add lib names
  cov=matrix(unlist(foreach(lib=all.counts.2alleles) %do% ((lib[,1]+lib[,2]))), ncol=length(all.counts.2alleles)) # coverage
  minfreqs[cov<mincov]=NA
  
  # all data: min allele freqs in base and hot
  if (is.null(polaRise))
  {
    data.table(data.frame(chr=x[,1], pos=x[,2], ref=x[,3], minallele, majallele, minfreqs))
  } else {

    if (is.list(polaRise))
    {
      afcs <- foreach(repl=polaRise, .combine=rbind) %do% {
        minfreqs[,1+repl][,2]-minfreqs[,1+repl][,1]
      }
      swap <- colMeans(afcs)<0

      riseallele <- ifelse(swap, majallele, minallele)
      fallallele <- ifelse(swap, minallele, majallele)

      cnames <- colnames(minfreqs)
      risefreqs <- foreach(lib=minfreqs, .combine=cbind) %do% {
        ifelse(swap, 1-lib, lib)
      }
      colnames(risefreqs) <- cnames

      data.table(data.frame(chr=x[,1], pos=x[,2], ref=x[,3], riseallele, fallallele, risefreqs))
    } else {
        stop("'ploaRise' has to be a list contaning number pairs, e.g. list(c(1,4),c(5,8),c(9,12))
                 with the number pairs describing the column numbers of the replicates between which
                 the allele frequency change (increase) is measured.")
    }

  }
  
  
  
  
  
}

#' @importFrom utils read.table
# check validity of the input parameters and return read in file
validity_sync_to_frequencies<-function(file, base.pops, header, mincov)
{
  # check parameter preence
  if (is.null(file)) { stop("Parameter 'file' missing with no default.") }
  if (is.null(base.pops)) { stop("Parameter 'base.pops' missing with no default.") }
  if (is.null(header)) { stop("Parameter 'header' missing with no default.") }
  if (mincov<10) { stop("Minimum coverage cannot be below 10.")
  } else if (mincov<15) { warning("Minimum coverages below 15 are not recommended due to unreliable frequency estimation.")}
  
  x=read.table(file,header=header)
  if (ncol(x)-3 != length(base.pops)) { stop("The number of libraries in the input sync file 'file' and the logical vector base.pops
                                             needs to be identical.")}
  return(x)
}




#' @import foreach
#' @importFrom stringi stri_split_fixed
#function which parses sync A:C:T:G:N:Indel format to matrix of ACTG
#d = sync file
#cols = columns to parse
getATCG <- function(d, cols)
{
  i <-NULL # declared in subsequent foreach
  foreach(i=cols) %do% 
  {
    cat(paste("getATCG of lib",colnames(d)[i],"\n"))
    matrix(as.numeric(unlist(stri_split_fixed(d[, i], pattern=":")))
           , nrow=nrow(d), byrow=TRUE)[, 1:4]
  }
}


#' @import data.table
#' @importFrom matrixStats rowRanks
#' @importFrom stats rnorm
# extends rowRanks function but ensures that with equal counts the same rank will not be given twice
# it adds to each count a small random value, so exact "counts" will never be the same
rowRanks2 <- function(v)
{
  vv <- v + (foreach(1:4, .combine=cbind) %do% rnorm(mean=0, sd=0.1, n=nrow(v)))
  rowRanks(vv)
}


#' @importFrom matrixStats rowRanks
#function which determines the two most common SNPs
#
top2 <- function(count.sum)
{
  topSNPs <- t(rowRanks2(count.sum))
  m3 <- which(topSNPs==3) %% 4 # row of the minor allele
  m3[m3==0] <- 4; if (sum(m3==0)>0) print(paste("m3",which(m3==0)))
  m4 <- which(topSNPs==4) %% 4 # row of major allele
  m4[m4==0] <- 4; if (sum(m4==0)>0) print(paste("m4",which(m3==0)))
  return(data.frame(Min=m3, Max=m4))
}
  
  
