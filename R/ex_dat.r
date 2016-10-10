
# ex_dat=sync_to_frequencies(file="../../../../haplotype_simul_number/MS03_haplo_reconstruct/03_simlulate_glm_2sites/sim58_2s.sync", base.pops=rep(c(TRUE,rep(FALSE,8)),times=5), header=FALSE)
# setnames(ex_dat, names(ex_dat)[7:ncol(ex_dat)], c("F0_R1","F10_R1","F20_R1","F30_R1","F40_R1","F50_R1","F60_R1","F70_R1","F80_R1",
#                                   "F0_R2","F10_R2","F20_R2","F30_R2","F40_R2","F50_R2","F60_R2","F70_R2","F80_R2",
#                                   "F0_R3","F10_R3","F20_R3","F30_R3","F40_R3","F50_R3","F60_R3","F70_R3","F80_R3",
#                                   "F0_R4","F10_R4","F20_R4","F30_R4","F40_R4","F50_R4","F60_R4","F70_R4","F80_R4",
#                                   "F0_R5","F10_R5","F20_R5","F30_R5","F40_R5","F50_R5","F60_R5","F70_R5","F80_R5"))
# ex_dat=ex_dat[,.(chr,pos,ref,minallele,majallele,basePops,
#           F0_R1,F20_R1,F40_R1,F60_R1,
#           F0_R2,F20_R2,F40_R2,F60_R2,
#           F0_R3,F20_R3,F40_R3,F60_R3,
#           F0_R4,F20_R4,F40_R4,F60_R4,
#           F0_R5,F20_R5,F40_R5,F60_R5)]
# # save(ex_dat, file="data/ex_dat.rda", compress = TRUE, compression_level=9)
# # data set slightly reduced as data for whole chr arm was too large for submission
# # to CRAN
# ex_dat <- ex_dat[6700000<=pos & pos<=16500000]
# save(ex_dat, file="data/ex_dat.rda", compress = TRUE, compression_level=9)


#' Example data set
#' 
#' Example data set including the simulated data set with
#' selection oparating on two different sites 1 Mb apart from 
#' each other, each unique to a single but different out
#' of 200 different founder haplotypes (simulated data 
#' corresponding to Fig. 1A in Franssen, Barton & Schloetterer 2016, 
#' \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}).
#' Data was simulated in an evolve and resequence setup
#' using 5 replicates and time points up to generation
#' F60 with sampled time-points every 20th generation. The
#' data as a data.table containing already allele frequencies
#' polarized for the minor allele in the founder populations
#' as can be obtained by the \code{\link{sync_to_frequencies}} function
#' from a sync formatted text file (sync format, see Kofler 
#' et al. 2011).
#' A basic workflow for haplotype analysis is presented 
#' below. 
#'
#' @name ex_dat
#' @docType  data
#' @usage data(ex_dat)
#' @format data.table containing frequency information of 
#' all samples polarized for the minor allele in the experimental
#' founder population
#' @references Franssen, Barton & Schloetterer 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @author Susanne U. Franssen
#' @seealso \code{\link{hbr}} 
#' \code{\link{sync_to_frequencies}} \code{\link{initialize_SNP_time_series}} \code{\link{SNP_time_series}}
#' \code{\link{reconstruct_hb}} \code{\link{plot.hbr}} \code{\link{summary.hbr}} \code{\link{plot_hbr_freq}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{{plot_cluster_trajectories}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}} 
#' \code{\link{inspect_window_dbScan}} \code{\link{markers}}\code{\link{number_hbr.hbr}}
#'
#' @examples
#' # The following workflow provides a general example on 
#' # haplotype analysis on the provided example data set
#' #
#' # example data was previously formated from a sync file with an added header using:
#' # ex_dat=sync_to_frequencies(file="ex_dat.sync", 
#' #                  base.pops=rep(c(TRUE,rep(FALSE,3)),times=5), header=TRUE)
#' # The file contains samples for F0, F20, F40 and F60, each for five replicates simualations.
#' 
#' # filter replicated time series data for informative SNPs
#' dat_filtered=initialize_SNP_time_series(chr=ex_dat$chr, pos=ex_dat$pos, 
#' base.freq=ex_dat$basePops, lib.freqs=ex_dat[,7:ncol(ex_dat), with=FALSE], 
#' pop.ident=rep(1:5,each=4), pop.generation=rep(c(0:3)*20,times = 5), use.libs=rep(TRUE,20))
#' 
#' # reconstruct haplotype-blocks
#' dat_reconst=reconstruct_hb(dat_filtered, chrom="2R")
#' 
#' # various ways of inspecting the results
#' #
#' ?plot.hbr
#' plot(dat_reconst, indicate_shared=TRUE, addPoints=TRUE)
#' #
#' #?summary.hbr
#' summary(dat_reconst)
#' #
#' ?plot_hbr_freq
#' par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
#' plot_hbr_freq(dat_reconst, hbr_id=2, replicate=1, timepoint=c(0,20,40,60), window=5)
#' plot_hbr_freq(dat_reconst, hbr_id=2, replicate=2, timepoint=c(0,20,40,60), window=5)
#' # Note: For the example parameter settings reconstructed haplotype-block hbr_id=2
#' # corresponds to the focal selected haplotype shown in Fig. 1A 
#' # (Franssen, Barton & Schloetterer 2016,
#' # Reconstruction of haplotype-blocks selected during experimental evolution, MBE).
#' #
#' ?map
#' map(dat_reconst)
#' #
#' ?rev_map
#' rev_map(dat_reconst)
#' #
#' ?plot_cluster_trajectories
#' plot_cluster_trajectories(dat_reconst, window=38)
#' #
#' ?plot_marker_trajectories
#' plot_marker_trajectories(dat_reconst, hbr_id=2)
#' #
#' ?inspect_window
#' inspect_window(dat_reconst, window=38)
#' #
#' ?inspect_window_PCA
#' inspect_window_PCA(dat_reconst, window=38)
#' #
#' ?inspect_window_avLink
#' inspect_window_avLink(dat_reconst, window=38)
#' #
#' ?inspect_window_dbScan
#' inspect_window_dbScan(dat_reconst, window=38, eps=1)
#' #
#' ?markers
#' markers(dat_reconst, hbr_id=2)
#' #
#' ?number_hbr
#' number_hbr(dat_reconst)
NULL









