


# setOldClass("data.table") tells the formal (S4) class system that you'd like to 
# treat an object with a class attribute data.table (S3 class) as an S4 class. 
#--> so it can be used as a slot in the following S4 classes
setOldClass("data.table")

#' An S4 class storing time series data
#' 
#' Time-series data is initialzed to be used for haplotype-block reconstruction.
#'
#' @aliases SNP_time_series
#' @author Susanne U. Franssen
#' @description An S4 class representing the input data for haplotype-block 
#' reconstruction
#' @details An S4 class storing the initialized input data for haplotype-block 
#' reconstruction with \code{reconstruct_hb}. Genome-wide time series data for 
#' multiple time points and replicates is stored along with replicate and time 
#' point information of the time series data. The stored data is filtered for
#' SNPs with a maximum minor frequency in the experimental starting population 
#' and a minimum frequency change \code{minfreqchange} that is required for 
#' \code{minrepl} many replicates. An object of the class 
#' \code{SNP_time_series} can only be created with the function 
#' \code{\link{initialize_SNP_time_series}}.
#' @import data.table
#' @slot col.info data.table with columns 'chr', 'pos', 'base.freq' and window. 
#' Each row corresponds to a SNP position that fullfills the filtering criteria.
#' @slot lib.freqs data.table with columns for the different libraries (time 
#' points and replicates) and rows for all SNP positions that fullfill the 
#' filtering criteria.
#' @slot pos.cM numeric vector corresponding to SNP positions in col.info with
#' genetic positions in cM.
#' @slot  pop.ident numeric vector specifying the identity of each library in 
#' terms of replicate ID
#' @slot  pop.generation numeric vector specifying the time point of the 
#' respective library
#' @slot use.libs logical vector specifying which libraries should be used for 
#' haplotype-block reconstruction
#' @slot  winsize numeric specifying the window size on which to perform 
#' the analysis
#' @slot win.scale character string specifying which genome-wide distance 
#' measure is used for window definition. Options are "bp" (base pairs) or 
#' "cM" (centi Morgan). cM distances can only be used if gentic positions are 
#' provided in 'pos.cM' (default="Mb").
#' @slot  min.minor.freq numeric specifying the minimum frequency of the 
#' minor allele (polarized in the experimental starting population) to be 
#' included in the analysis (default=0).
#' @slot  max.minor.freq numeric specifying the maximum frequency of the 
#' minor allele (polarized in the experimental starting population) to be 
#' included in the analysis (default=3/113).
#' @slot min.lib.frac minimum fraction of non-NA values for a SNP across 
#' libraries (only using libraries specified in \code{use.libs}) (default=0.75).
#' @slot  minfreqchange numeric specifying the minimum frequency change 
#' required in 'minrepl' replicates 
#' required to
#' include the SNP in the analysis
#' @slot  minrepl numeric specifying the number of replicates, in which 
#' the 'minfreqchange' is required 
#' to include the SNP in the analysis
#' @references Franssen, Barton & Schloetterer 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @seealso \code{\link{ex_dat}} \code{\link{initialize_SNP_time_series}} \code{\link{reconstruct_hb}} 
#' 
setClass("SNP_time_series",
         representation(col.info="data.table",
                        lib.freqs="data.table",
                        pos.cM="numeric",
                        
                        pop.ident="numeric",
                        pop.generation="numeric",
                        use.libs="logical",
                        
                        winsize="numeric",
                        win.scale="character",
                        
                        min.minor.freq="numeric",
                        max.minor.freq="numeric",
                        min.lib.frac="numeric",
                        
                        minfreqchange="numeric",
                        minrepl="numeric"))





#' An S4 class storing results from haplotype-block reconstruction
#'
#' @aliases hbr
#' @author Susanne U. Franssen
#' @description An S4 class containg input parameters and results for reconstructed 
#' haplotype-blocks
#' @details An S4 class containing a SNP_time_series data object, the haplotype-blocks were
#' constructed on, parameters that were used for haplotype-block reconstruction and results
#' of haplotype block reconstruction. An hbr object can only be created with 
#' \code{\link{reconstruct_hb}}.
#' @import data.table
#' @slot dat A SNP_time_series object containg the time series data and threshold used for 
#' marker filtering (refer to \code{SNP_time_series} for more details)
#' @slot chromosome The chromosome for which haplotype-blocks were calculated.
#' @slot min.cl.size Numeric specifying the minimum number of correlated markers in a 
#' window for haplotype-block reconstruction
#' @slot min.cl.cor Numeric specifying the correlation between marker SNPs required 
#' using average linkage clustering for markers to be assembled in one cluster
#' @slot  min.inter Numeric specifying the minimum number of markers in two clusters 
#' of overlapping windows required to be identical for cluster elongation across windows 
#' to build haplotype-blocks
#' @slot single.win Boolean specifying that are supported by only a cluster in one window are 
#' also included in the markers of reconstructed blocks. If FALSE only hbrs are returned that span
#' at least two windows and only markers being present in the intersection between overlapping
#' windows are included.
#' @slot transf Boolean indicating if time series data was sqrt transformed
#' prior to clustering.
#' @slot arcsine Boolean indicating if time series data was arcsine transformed
#' prior to clustering.
#' @slot scaleSNP Boolean indicating if time series data was scaled (mean=0, var=1)
#' for each SNP prior to clustering.
#' @slot pos.cor Boolean indicating if negative correletions between time series between
#' two SNPs were set to zero prior to clustering.
#' @slot cl.chr A list containing results from clustering of low frequency markers in the 
#' experimental starting population for each window. Each list element contains the results 
#' for a window.If clusters for the respective window were identified the respective list
#' numeric vectors with the marker positions for each cluster. If no clusters were identified 
#' the list is empty. The clusters identified for overlapping windows are the basis for 
#' cluster elongation to haplotype-blocks as present in \code{cl.long.m}.
#' @slot cl.long.m A list containing the results of haplotype-bock reconstruction. Each 
#' list element corresponds to a reconstructed #' haplotype-block and contains a numeric 
#' vector with all the SNP positions. The minor allele in the experimnental starting
#' population at those positions represent marker alleles for the respective hapltotype-block. 
#' Chromosome-wide haplotype-blocks can be visualized with \code{plot}.
#' @references Franssen, Barton & Schloetterer 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @seealso \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'
setClass("hbr",
         representation(dat="SNP_time_series",
                        
                        chromosome="character",

                        min.cl.size="numeric",
                        min.cl.cor="numeric",
                        
                        min.inter="numeric",
                        single.win="logical",
                        
                        transf="logical",
                        arcsine="logical",
                        scaleSNP="logical",
                        pos.cor="logical",
                        
                        cl.chr="list",
                        cl.long.m="list"))



