

#' Display of a hbr object
#'
#' This function displays the summarized features of a hbr object
#'
#' @aliases show show.hbr
#' @author Susanne U. Franssen
#' @details This function displays the summarized features of a hbr object
#' @export
#' @importFrom methods show
#' @param object object of the class hbr data
#' @references Franssen et al. 2016, \href{http://mbe.oxfordjournals.org/content/early/2016/10/03/molbev.msw210.abstract}{Reconstruction of haplotype-blocks
#' selected during experimental evolution}, \href{http://mbe.oxfordjournals.org/}{MBE}
#' @seealso \code{\link{summary.hbr}} \code{\link{ex_dat}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#' 
setMethod("show","hbr",
          function(object) {
            cat("Object of class",class(object),"\n\n")
            
            cat("Information on used time series data:\n")
            cat("Number of populations included:",length(unique(object@dat@pop.ident)),"\n")
            cat("Number of libraries included:",length(object@dat@pop.ident),"\n")
            cat("Number of libraries used for filtering:",sum(object@dat@use.libs),"\n")
            cat("Data has been filtered for SNPs with:\n"
                ,"\t allele frequencies >=",object@dat@min.minor.freq," in the founder population\n"
                ,"\t allele frequencies <=",object@dat@max.minor.freq," in the founder population\n"
                ,"\t a frequency increase of",object@dat@minfreqchange
                ," in at least ",object@dat@minrepl,"replicates\n")
            cat("Number of SNPs after filtering:",nrow(object@dat@col.info),"\n")
            cat("Window definition:\n")
            cat("\t distance measure used:",object@dat@win.scale,"\n")
            cat("\t applied windowsize:",object@dat@winsize,"\n")
            
            cat("Parameters used for haplotype-block reconstruction:\n")
            cat("Chromosome",object@chromosome,"\n")
            cat("Parameters used for each sliding window:\n"
                ,"\t a minimum cluster size of ", object@min.cl.size,"\n"
                ,"\t an average correlation within a cluster of ", object@min.cl.cor,"\n")
            cat("Parameters for cluster elongation across windows:\n"
                ,"\t a minimum number of ",object@min.inter,"intersecting markers\n")
            cat("single.win:",object@single.win,"\n")
            cat("transf:",object@transf,"\n")
            cat("arcsine:",object@arcsine,"\n")
            cat("scaleSNP:",object@scaleSNP,"\n")
            cat("pos.cor:",object@pos.cor,"\n\n")
            
            
            cat("Summary of results:\n")
            cat("Window summary:\n"
                ,"\t for ",sum(unlist(lapply(object@cl.chr,length))!=0)
                ," windows clusters were identified \n")
            #cat(", with\n\t window_id #_cluster \n")
            #print_wincl_hbr_info(object@cl.chr)
            cat("Haplotype-block (hbr) summary:\n"
                ,"\t",length(object@cl.long.m)," haplotype-blocks (hbr) were reconstructed\n")
            #cat(", with\n\t hbr_id #_marker \n")
            #print_wincl_hbr_info(object@cl.long.m)
            cat("\t hbrs (and markers) where included if they\n"
                ,"\t\t are only supported by a single window:"
                ,object@single.win)
          })





#' Method to summarize information of reconstructed haplotype-blocks
#'
#' The method summarizes information of the reconstructed haplotype-blocks for a chromosome.
#'
#' @aliases summary.hbr
#' @author Susanne U. Franssen
#' @details The method operates on \code{\link{hbr}} objects and summarizes information of 
#' reconstructed haplotype-blocks
#' @export
#' @param object object of class \code{\link{hbr}} containing the results of reconstructed haplotype
#' blocks to be summarized.
#' @param min.marker numeric specifying the minimum number of markers a haplotype-block to be 
#' reported. NOTE: IDs of haplotype-blocks in the provided summary are identical to the 
#' IDs of the previously identified blocks. 
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#' @return data.table with row entries for each reconstructed haplotype-block and columns for:
#' chr (chromosome), id (numbering of selected blocks), n.marker (number of markers in the 
#' reconstructed block), spos, epos (starting and end position of the block), len.pos (length
#' of the block in bp), MperMb (number of markers per Mb), swin, ewin (starting and end window of the block)
#'
setMethod("summary","hbr",
          function(object, min.marker=1) {
            
            mil=1000000
            len=length(object@cl.long.m)
            chr=rep(object@chromosome,len)
            id=1:len
            n.marker=unlist(lapply(object@cl.long.m,length))  # number of markers
            spos=unlist(lapply(object@cl.long.m,min))         # start pos
            epos=unlist(lapply(object@cl.long.m,max))         # end pos
            swin=lapply( map(object), FUN = function(x) x[1]) # start window
            ewin=lapply( map(object), FUN = function(x) x[length(x)]) # end window
            len.pos=epos-spos+1                               # lenght in bp
            MperMb=round(n.marker/len.pos*mil,1)              # marker per Mb
            
            res=data.table(chr,id,n.marker,spos,epos,len.pos,MperMb)
            res$id=as.integer(id)
            res$n.marker=as.integer(n.marker)
            res$spos=as.integer(spos)
            res$epos=as.integer(epos)
            res$swin=as.integer(swin)
            res$ewin=as.integer(ewin)
            res$len.pos=as.integer(len.pos)
            res$MperMb=as.double(MperMb)
            if (min.marker>1)
            {
              res=res[n.marker>=min.marker]
              # update id
              #res$id=1:nrow(res)
            }

            return(res)
          })


#' Method to visualize reconstructed haplotype-blocks
#'
#' The method visualizes reconstructed haplotype-blocks for a chromosome.
#'
#' @aliases plot.hbr
#' @author Susanne U. Franssen
#' @details The method operates on \code{\link{hbr}} objects and visualizes location of 
#' reconstructed haplotype-blocks with respect to its genomic position.
#' @export
#' @importFrom graphics segments abline
#' @importFrom utils combn
#' @param x object of class \code{hbr} containing the results of reconstructed haplotype
#' blocks for visualization.
#' @param min.marker numeric specifying the minimum number of markers a haplotype-block 
#' should contain in order to be visualized.
#' @param xlab Label of the x-axis with the default value 'Genomic position [MB]'.
#' @param ylab Label of the y-axis with the default value 'Reconstructe haplotype-block'.
#' @param main Plot title (default: "Chromosome XX").
#' @param col Color of the lines representing the haplotype-blocks (default: "black").
#' @param lwd Line width of the lines representing the haplotype-blocks (default: 4).
#' @param hline Distance between horizontal lines plotted for orientation (default: 10).
#' @param ... arguments of the generic plot method.
#' @param indicate_shared logical value specifying if "spurious" markers that are 
#' identical between pairs of haplotype-blocks should be indicated.
#' This function is usefull for inspecting results and deciding whether all
#' identified blocks are all independent or maybe reconstruction parameters should 
#' be changed. 
#' @param addPoints logical value indicating if for each reconstructed block markers
#' should additionally be indicated.
#' @param hbr_plot boolean vector of length the number of reconstructed bocks indicating which
#' ones should be plotted
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}}
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}}\code{\link{number_hbr}}
#'
setMethod("plot","hbr",
          function(x, min.marker=1, xlab="Genomic position [Mb]", 
                   ylab="Reconstructed haplotype-block", main="default",
                   col="black", lwd=4, hline=10, 
                   indicate_shared=F, addPoints=F,
                   hbr_plot=NULL, ...)
          {
            if (is.null(hbr_plot))
            {
              hbr_plot=rep(T,length(x@cl.long.m))
            } else if (length(hbr_plot) != length(x@cl.long.m))
            {
              stop(paste("hbr_plot must be a boolean vector of length",length(x@cl.long.m)))
            }
            hbr.i=unlist(lapply(x@cl.long.m,length))>=min.marker # boolean for hbrs with sufficient markers
            hbr.i=(hbr_plot+hbr.i==2) # only hbrs with enough min markers and were specified to plot
            hbrs=x@cl.long.m[hbr.i] # hbrs of sufficient # of markers
            
            mil=1000000
            min.pos=min(unlist(hbrs))/mil
            max.pos=max(unlist(hbrs))/mil
            
            if(main=="default"){main=paste("Chromosome",x@chromosome)}
            if(length(col)<length(hbrs)){col=rep(col,length(hbrs))}
            plot(x=c(min.pos,max.pos),y=c(0.5,length(hbrs)+0.5), col="white", xlab=xlab, ylab=ylab, main=main, ...)
            abline(h=(0:round(length(hbrs)/hline))*hline, lwd=0.7, lty=2)
            
            for (i in 1:length(hbrs))
            {
              if (addPoints) points(hbrs[[i]]/mil, rep(i,length(hbrs[[i]])), col="red")
              segments(min(hbrs[[i]])/mil,i,max(hbrs[[i]])/mil,i,lwd=lwd, col=col[i])
            }
            
            if (indicate_shared & length(hbrs)>1)
            {
              # all pairwise comparisons of marker sets of hbrs
              hbrNames.comb=combn( 1:length(hbrs) , 2 , simplify = FALSE )
              hbrSets.comb=combn( hbrs , 2 , simplify = FALSE )
              hbrinter.comb=lapply(hbrSets.comb , function(x) intersect(x[[1]], x[[2]]) )
              for (comp in 1:length(hbrinter.comb))
              {
                #draw connecting lines whenever markers are shared between two hbrs
                if (length(hbrinter.comb[[comp]]) >0)
                {
                  segments(x0=hbrinter.comb[[comp]]/mil, y0=hbrNames.comb[[comp]][1], 
                           x1=hbrinter.comb[[comp]]/mil, y1=hbrNames.comb[[comp]][2])
                }
              }
            }
          })


#' Method to visualize the trajectories for all identified clusters in a window
#'
#' @aliases plot_cluster_trajectories plot_cluster_trajectories.hbr
#' @author Susanne U. Franssen
#' @details The method operates on \code{\link{hbr}} objects and plots the trajectories
#' of all clusters and replicates in the specified window. Cluster trajectories are
#' visualized by trajectories of all markers contained in the respective clusters and
#' median trajectories for a cluster are visualized in a different color. For each
#' cluster it is indicated to which haplotype-block it was later assigned.
#' NOTE: It is possible that to clusters are assigned to the identical
#' haplotype-block. This can happen when a cluster in an overlapping window shares 
#' markers with two clusters of the "focal" window.
#' @export
#' @importFrom grDevices rainbow
#' @importFrom stats median
#' @param object object of class \code{hbr}
#' @param window numeric specifying the number of the window, for which trajectories 
#' should be plotted.
#' @param ylim numeric vector with two elements specifying the limits of the y-axis.
#' @param p.median boolean specifying if the median trajectory for each cluster should 
#' be added to the plot (default: TRUE).
#' @param tp.cor Boolean indicating if only the time points used for calaculating
#' correlations \code{use.libs} are shown (tp.cor=T) or all time points present in the
#' data set are shown (tp.cor=F).
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'

setMethod("plot_cluster_trajectories", "hbr",
          function(object, window, ylim=c(0,1), p.median=T, tp.cor=T)
          {
            if (!window %in% 1:length(object@cl.chr))
            {
              writeLines(paste("This window does not exist. Windows range from 1 to",length(object@cl.chr),"."))
            } else if (length(object@cl.chr[[window]])==0)
            {
              writeLines("No clusters were called for this window.")
            } else # print trajectories for all clusters identified for this window
            {
              # for adding hbr info
              revMap=rev_map(object)
              hbrs=unlist(revMap[names(revMap) == paste0("w",window)]) 
              
              n.cluster=length(object@cl.chr[[window]]) # number of clusters identified for this window
              n.repl=length(unique(object@dat@pop.ident)) # number of different replicates / populations provided
              xmin=min(object@dat@pop.generation)
              xmax=max(object@dat@pop.generation)
              cols=rainbow(n.cluster) # colors for median lines in different clusters
              
              par(mfrow=c(n.cluster,n.repl),mar=c(0,0,0,0),oma=c(4,4,3,1))
              for (i in 1: n.cluster)
              {
                SNP.pos=object@cl.chr[[window]][[i]]
                # test which hbr this cluster belongs to
                hbr=c()
                for (xx in hbrs)
                {
                  if (length(intersect(object@cl.long.m[[xx]], SNP.pos)) 
                      >= object@min.inter)#(0.2)*length(SNP.pos))
                  { hbr=c(hbr,xx)}
                }
                
                pos <- NULL # silly as pos is defined in object@dat@col.info$pos
                if (tp.cor) # only time points used for the correlations will be shown
                {
                  lf=object@dat@lib.freqs[object@dat@col.info$pos %in% SNP.pos,object@dat@use.libs,with=F] #library frequencies
                } else {
                  lf=object@dat@lib.freqs[object@dat@col.info$pos %in% SNP.pos,] #library frequencies
                }
                bf=object@dat@col.info[pos %in% SNP.pos]$base.freq #composite base frequencies 
                if(i==n.cluster){xaxt=NULL}else{xaxt="n"}
                for (repl in 1:n.repl) # plot for each replicate population
                {
                  if (tp.cor) # only time points used for the correlations will be shown
                  {
                    gens=object@dat@pop.generation[object@dat@use.libs][repl==object@dat@pop.ident[object@dat@use.libs]]
                    lf.repl=lf[,repl==object@dat@pop.ident[object@dat@use.libs],with=F][,order(gens),with=F]
                  } else {
                    gens=object@dat@pop.generation[repl==object@dat@pop.ident]
                    lf.repl=lf[,repl==object@dat@pop.ident,with=F][,order(gens),with=F]
                  }
                  
                  if (sum(gens==0) == 0) # add composite base to trajectory
                  {
                    lf.repl=cbind(bf,lf.repl)
                    gens=c(0,gens)
                  }
                  if(repl==1){yaxt=NULL}else{yaxt="n"}
                  plot(1,col="white", xlim=c(xmin,xmax),ylim=ylim,xaxt=xaxt,yaxt=yaxt)
                  text(x=xmax*0.2,y=ylim[2]*0.9, paste("hbr",paste(hbr,collapse=","))) # hbr info
                  
                  for(line in 1:nrow(lf.repl))
                  {
                    lines(gens,lf.repl[line,])
                  }
                  if(p.median){lines(gens,apply(lf.repl,2,median),col=cols[i],lwd=1.5)}
                }
              }
              mtext("Generation",1,outer = T,line = 3)
              mtext("Marker frequency",2,outer = T,line = 3)
              mtext(paste("Window",window," Replicates",paste(sort(unique(object@dat@pop.ident)),collapse=", ")),3,outer = T,line = 1)
              
            }
            
          })


#' Method to visualize the trajectories for all markers in a haplotype-block
#'
#' @aliases plot_marker_trajectories plot_marker_trajectories.hbr
#' @author Susanne U. Franssen
#' @details The method operates on \code{\link{hbr}} objects and plots the trajectories
#' of all markers in a haplotype-block in all replicates.
#' Note: As blocks can span a wide range along the genome, it is not expected that 
#' trajectories within a block stay very similar. Changing location along the genome 
#' can be indicated with loc.col=T.
#' @export
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics legend lines mtext par points text
#' @param object object of class \code{hbr}
#' @param hbr_id the id of the haplotype-block, for which the marker trajectories 
#' should be plotted. 
#' @param ylim numeric vector with two elements specifying the limits of the y-axis.
#' @param loc.col boolean indicating if trajectories should be coloured ranging from
#' red over blue to yellow according to their location on the chromosome (default: TRUE).
#' @param tp.cor Boolean indicating if only the time points used for calaculating
#' correlations \code{use.libs} are shown (tp.cor=T) or all time points present in the
#' data set are shown (tp.cor=F).
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'

setMethod("plot_marker_trajectories", "hbr",
          function(object, hbr_id, ylim=c(0,1), loc.col=T, tp.cor=T)
          {
            marker=markers(object, hbr_id)
            n.repl=length(unique(object@dat@pop.ident)) # number of different replicates / populations provided
            xmin=min(object@dat@pop.generation)
            xmax=max(object@dat@pop.generation)
            if (loc.col)
            {
              colfunc <- colorRampPalette(c("red", "blue","yellow"))
              cols=colfunc(length(marker))
            } else {cols=rep("black",length(marker))}
              
              
            par(mfrow=c(1,n.repl),mar=c(0,0,0,0),oma=c(4,4,3,1))

            SNP.pos=marker

            pos <- NULL # silly as pos is defined in object@dat@col.info$pos
            if (tp.cor) # only time points used for the correlations will be shown
            {
              lf=object@dat@lib.freqs[object@dat@col.info$pos %in% SNP.pos,object@dat@use.libs,with=F] #library frequencies
            } else {
              lf=object@dat@lib.freqs[object@dat@col.info$pos %in% SNP.pos,] #library frequencies
            }
            bf=object@dat@col.info[pos %in% SNP.pos]$base.freq #composite base frequencies 
            for (repl in 1:n.repl) # plot for each replicate population
            {
                 if (tp.cor) # only time points used for the correlations will be shown
                 {
                   gens=object@dat@pop.generation[object@dat@use.libs][repl==object@dat@pop.ident[object@dat@use.libs]]
                   lf.repl=lf[,repl==object@dat@pop.ident[object@dat@use.libs],with=F][,order(gens),with=F]
                 } else {
                   gens=object@dat@pop.generation[repl==object@dat@pop.ident]
                   lf.repl=lf[,repl==object@dat@pop.ident,with=F][,order(gens),with=F]
                 }
                 
                 if (sum(gens==0) == 0) # add composite base to trajectory
                 {
                   lf.repl=cbind(bf,lf.repl)
                   gens=c(0,gens)
                 }
                 if(repl==1){yaxt=NULL}else{yaxt="n"}
                xaxt=NULL
                 plot(1,col="white", xlim=c(xmin,xmax),ylim=ylim,xaxt=xaxt,yaxt=yaxt)
                 
                 for(line in 1:nrow(lf.repl))
                 {
                   lines(gens,lf.repl[line,], col=cols[line])
                 }
            }
            mtext("Generation",1,outer = T,line = 3)
            mtext("Marker frequency",2,outer = T,line = 3)
            
          })


#' Map from reconstructed haplotype-blocks to windows
#' 
#' A method returning information which windows are contain in which reconstructed
#' haplotype-block.
#' 
#' @aliases map map.hbr
#' @author Susanne U. Franssen
#' @details The method operates on \code{\link{hbr}} objects and returns a list summarizing
#' the mapping information from reconstructed haplotype-blocks to contained windows.
#' NOTE: A window is only listed for a haplotype-block when at least \code{min.inter}
#' many markers of the block are in that window.
#' @export
#' @param object object of class \code{\link{hbr}} containing the results of reconstructed 
#' haplotype-blocks
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#' @return list with elements corresponding to all reconstructed haplotype-blocks (id for all
#' blocks without filtering based on the number of markers). Each element contains a vector of
#' window indices that are contained in the block
#'
setMethod("map", "hbr",
          function(object)
          {
            hb <- NULL # initialized in the subsequent foreach
            win.list = foreach (hb=object@cl.long.m) %do%
            {
              wins=c()
              for (win in 1:length(object@cl.chr))
              {
                # intersection between blocks in that marker and in any cluster of the respective 
                # window
                inter=intersect(hb,unlist(object@cl.chr[[win]]))
                if (length(inter)>=object@min.inter) { wins=c(wins,win) }
              }
              #paste(wins,collapse = ",")
              wins
            }
            
            # each list element corresponds to an hbr and contains a vector of all windows it spans
            win.list  
          })


#' Map from windows to reconstructed haplotype-blocks
#' 
#' A method returning information which reconstructed haplotype-blocks are present in 
#' which windows. 
#' 
#' @aliases rev_map rev_map.hbr
#' @author Susanne U. Franssen
#' @details The method operates on \code{\link{hbr}} objects and returns a list summarizing
#' the mapping information from windows to reconstructed haplotype-blocks.
#' @export
#' @param object object of class \code{\link{hbr}} containing the results of reconstructed 
#' haplotype-blocks
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#' @return list with elements corresponding to all windows that contain a reconstructed
#' haplotype-block segment. Each element contains a vector of haplotype-block indices
#' that overlap the respective window.
#'
setMethod("rev_map", "hbr",
          function(object)
          {
            reverse_map=list()
            xx=map(object)
            for (i in 1:length(xx)) # over all hbr
            {
              for (j in 1:length(xx[[i]])) # over all win in each hbr
              {
                if (!paste0("w",xx[[i]][j]) %in% names(reverse_map))
                {
                  reverse_map[[length(reverse_map)+1]]=c(i)
                  names(reverse_map)[length(reverse_map)]=paste0("w",xx[[i]][j])
                } else
                {
                  ii=(1:length(reverse_map))[names(reverse_map)==paste0("w",xx[[i]][j])]
                  reverse_map[[ii]]=c(reverse_map[[ii]],i)
                }
              }
            }
            reverse_map
            
          })



#' Markers of the specified reconstructed haplotype-block
#' 
#' Returns the genomic positions of the markers for the specified reconstructed haplotype-block
#' 
#' @aliases markers markers.hbr
#' @author Susanne U. Franssen
#' @export
#' @param object object of class \code{\link{hbr}} containing the results of reconstructed 
#' haplotype-blocks
#' @param hbr_id index of the haplotype-block of interest
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#' @return numeric vector with marker positions of the specified reconstructed haplotype-block
#'
setMethod("markers", "hbr",
          function(object, hbr_id)
          {
            if (hbr_id > length(object@cl.long.m)) stop(paste("Please provide an 'hbr_id' <=",length(object@cl.long.m),"."))
            object@cl.long.m[[hbr_id]]
          })

#' The number of reconstructed haplotype-block
#' 
#' Returns the number of reconstructed haplotype-blocks
#' 
#' @aliases number_hbr number_hbr.hbr
#' @author Susanne U. Franssen
#' @export
#' @param object object of class \code{\link{hbr}} containing the results of reconstructed 
#' haplotype-blocks
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}}
#' @return returns the number of reconstructed haplotype-blocks
#'
setMethod("number_hbr", "hbr",
          function(object)
          {
            length(object@cl.long.m)
          })



#' Plots frequencies of a reconstructed haplotype-block along the chromosome
#' 
#' @aliases plot_hbr_freq plot_hbr-freq.hbr
#' @author Susanne U. Franssen
#' @details The plotting method operates on \code{\link{hbr}} objects and returns a plot containing
#' frequencies of a respective haplotype-block along the chromosome for specified time points and 
#' replicates.
#' @export
#' @importFrom grDevices rainbow
#' @importFrom zoo rollapply
#' @param object object of class \code{\link{hbr}} containing the results of reconstructed 
#' haplotype-blocks
#' @param hbr_id id (integer) of the haplotype-block to be plotted
#' @param replicate numeric vector of integers specifying the replicates for which results
#' should be plotted
#' @param timepoint numeric vector of time points for which the results shouls be plotted
#' @param window window size over which frequencies are averaged (results are plotted for windows
#' overlapping by windows/2)
#' @param cols vector of colors with the length of the specified libraries
#' @param add logical specifying if a new plot should be created or haplotype-block frequencies 
#' are added to an existing plot.
#' @param sumstat Summary statistics used for the y-values in a window. Either 
#' specify "mean" (the default) or "median".
#' @param cex scaling of the point size in the output plot
#' @param xlab x-label in the output plot
#' @param ylab y-label in the output plot
#' @param xlim vector of the limits on the x-axis in the output plot
#' @param ylim vector of the limits on the y-axis in the output plot
#' @param pch option to specify symbols to use when plotting points in the output plot
#' @param lwd line width in the output plot
#' @param xaxt print "n" if x-axis scale should not be printed
#' @param yaxt print "n" if y-axis scale should not be printed
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} 
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'
setMethod("plot_hbr_freq", "hbr",
          function(object, hbr_id=1, replicate, timepoint, window=1, cols=NULL, add=F
                   , sumstat="mean", cex=0.7
                   , xlab="Genomic position [Mb]", ylab="Marker frequency", xlim=NULL
                   , ylim=c(0,1), pch=20, lwd=2
                   , xaxt=NULL, yaxt=NULL)
          {
            validity_plot_hbr_freq(object, hbr_id, replicate, timepoint, window
                                   , cols, add, cex
                                   , xlab, ylab, xlim, ylim)
            
            mil=1000000
            n.libs=length(object@dat@pop.ident)

            # get markers for the respective block
            # index of hbr markers in the time series data
            i.hbr.m=object@dat@col.info$pos %in% object@cl.long.m[[hbr_id]]
            col.info.hbr=object@dat@col.info[i.hbr.m,]
            lib.freqs.hbr=object@dat@lib.freqs[i.hbr.m,]
            if (nrow(col.info.hbr)<(1.5*window) & window>1){
              stop(paste("There are only",nrow(col.info.hbr),"markers, please provide a smaller 
                         window size."))
            }
            
            # get markers for the respective libraries
            libs.rep=(1:n.libs)[object@dat@pop.ident %in% replicate]
            libs.tp=(1:n.libs)[object@dat@pop.generation %in% timepoint]
            libs=intersect(libs.rep,libs.tp)
            if (length(libs)<1){
              stop("There are no libraries with the specified replicate AND timepoint properties.
                   Please provide valid values for replicate and timepoints.")}
            
            # color handling
            if (is.null(cols)){
              cols=rainbow(length(libs))
            } else {
              if (length(cols)!=length(libs)) {
                print(libs)
                stop(paste("Please provide a vector of",length(libs),"colors for the libraries",
                           paste(colnames(lib.freqs.hbr)[libs],collapse = ", ")))
              }
            }
            cat(paste("Block frequencies are provided for libraries:\n"
                        ,paste(colnames(lib.freqs.hbr)[libs],collapse = ", "),"\n"))
            
            # window smoothing
            if (window>1) {
              x=rollapply(col.info.hbr$pos[!is.na(lib.freqs.hbr[[libs[1]]])]
                          ,width=window,FUN=mean,by=1)/mil#window/2)/mil
              if (sumstat=="median"){
                y=rollapply(lib.freqs.hbr[[libs[1]]][!is.na(lib.freqs.hbr[[libs[1]]])]
                            ,width=window,FUN=median,by=1)#window/2)
              } else {
                y=rollapply(lib.freqs.hbr[[libs[1]]][!is.na(lib.freqs.hbr[[libs[1]]])]
                            ,width=window,FUN=mean,by=1)#window/2)
              }
              
            } else {
              x=col.info.hbr$pos/mil
              y=lib.freqs.hbr[[libs[1]]]
            }
            
            # plotting parameters
            if (is.null(xlim)){
              xlim=c(min(x),max(x))
            } else {xlim=xlim/mil}
            
            # plotting
            if (add){
              points(x,y,pch=pch,col=cols[1], cex=cex)
            } else {
              plot(x,y,pch=pch,col=cols[1], xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, cex=cex, xaxt=xaxt, yaxt=yaxt)
              abline(h=(0:5)*0.2,lty=2,lwd=0.5)
            }
            if (length(libs)>1)
            {
              for (i in 2:length(libs))
              {
                if (window>1){
                  x=rollapply(col.info.hbr$pos[!is.na(lib.freqs.hbr[[libs[i]]])]
                              ,width=window,FUN=mean,by=1)/mil#window/2)/mil
                  if (sumstat=="median"){
                    y=rollapply(lib.freqs.hbr[[libs[i]]][!is.na(lib.freqs.hbr[[libs[i]]])]
                                ,width=window,FUN=median,by=1)#window/2)
                  } else {
                    y=rollapply(lib.freqs.hbr[[libs[i]]][!is.na(lib.freqs.hbr[[libs[i]]])]
                                ,width=window,FUN=mean,by=1)#window/2)
                  }
                  
                } else {
                  y=lib.freqs.hbr[[libs[i]]]
                }
                points(x,y,pch=pch,col=cols[i],cex=cex)
              }
            }
          })


#' Haplotype-block marker inspection for a window
#' 
#' Plots correlations and groupings of haplotype-markers for a given window
#' 
#' @aliases inspect_window inspect_window.hbr
#' @author Susanne U. Franssen
#' @details The plotting method operates on \code{\link{hbr}} objects. For a 
#' specified window the correlation matrix for identified haplotype-block markers
#' in this window is visualized also indicating haplotype-block identity.
#' @export
#' @importFrom gplots heatmap.2
#' @importFrom grDevices colorRampPalette rainbow
#' @param object object of class \code{\link{hbr}} containing the results of 
#' reconstructed haplotype-blocks
#' @param window number of the window, which should be inspected
#' @param colCluster boolean value indicating whether columns in the output
#' correlation plot should be clustered or the original positions of the 
#' markers reflecting genomic positions should be kept.
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'
setMethod("inspect_window", "hbr",
          function(object, window, colCluster=T)
          {
            # extract data matrix from the respective window
            choose=( object@dat@col.info$win==window | object@dat@col.info$win==window+1) # relates to the data matrix in object
            sum(choose)
            ttpos=object@dat@col.info$pos[choose]
            #
            # data matrix in object@dat may contain some SNPs in between that were not included in any large cluster
            choose2= (object@dat@col.info$pos >= min(ttpos) &  object@dat@col.info$pos <= max(ttpos))# relates to the data matrix in object@dat
            aapos=object@dat@col.info$pos[choose2]
            #aadat=object@dat@lib.freqs[choose2,]
            aadat=object@dat@lib.freqs[choose2,object@dat@use.libs,with=F]
            rownames(aadat)<-aapos
            
            #
            if (object@transf) {aadat=sqrt(aadat)}#; print("trans")}
            if (object@scaleSNP) {aadat=data.table(t(scale(t(aadat))))}#; print("scale")}
            
            aadat_t=t(aadat); colnames(aadat_t)=aapos
            dwin.cor=cor(aadat_t,use="pairwise.complete.obs")
            if (object@pos.cor) {dwin.cor[dwin.cor<0]=0}
            
            # color according to performed clustering
            colV=rep("black", length(aapos))
            # block ids for window
            hbs=rev_map(object)[names(rev_map(object))==paste0("w",window)][[1]]
            cc=rainbow(length(hbs))
            for (i in hbs)
            {
              colV[ aapos %in% markers(object,i) ] <- cc[which(hbs==i)]
            }
            
            #require(gplots)
            colfunc <- colorRampPalette(c("yellow", "darkgreen"))
            
            # usage of average linkage clustering for heatmap
            hclust2 <- function(x, method="average")
              hclust(x, method=method)
            
            # heatmap colors
            # to keep the same colors across plots regardeles of the range of correlation values
            if(round(min(dwin.cor)/ (1/30)) == 30)
            {
              colH=colfunc(30)[29:30]
            } else {
              if (round(min(dwin.cor)<0)) {zz=0} else {zz=round(min(dwin.cor)/ (1/30))}
              colH=colfunc(30)[zz:30]
            }
            
            if (colCluster)
            {
              heatmap.2(dwin.cor, ColSideColors=colV, RowSideColors=colV, scale="none",
                        key.xlab="Correlation", hclustfun=hclust2, 
                        main=paste("Window",window, " hbr",paste(hbs,collapse = ",")), 
                        trace="none", key.title="", 
                        col=colH, dendrogram ="both")
            } else {
              heatmap.2(dwin.cor, ColSideColors=colV, RowSideColors=colV, scale="none", 
                        key.xlab="Correlation", hclustfun=hclust2,
                        main=paste("Window",window, " hbr",paste(hbs,collapse = ",")), 
                        trace="none", key.title="", 
                        col=colH, dendrogram ="row", Colv=F)
            }
          })





#' PCA of haplotype-block marker for a window
#' 
#' Plots all filtered SNPs for a given window in terms of their the 
#' first two principal components of the respective time series
#' data.
#' 
#' @aliases inspect_window_PCA inspect_window_PCA.hbr
#' @author Susanne U. Franssen
#' @details Performs principal component analysis (PCA) on the time series data 
#' of the filtered SNPs and the selected time points (filtered and selected 
#' through \code{use.libs} in \code{\link{initialize_SNP_time_series}})
#' in a given window. SNPs are shown in terms of 
#' the first two principal components and colored according to the
#' haplotype-block they were assigned to. Black colored points were 
#' not assigned to any cluster.
#' @export
#' @importFrom grDevices rainbow
#' @importFrom stats prcomp
#' @param object object of class \code{\link{hbr}} containing the results of 
#' reconstructed haplotype-blocks
#' @param window number of the window, which should be inspected
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_avLink}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'
setMethod("inspect_window_PCA", "hbr", 
          function(object, window)
          {
            # extract data matrix from the respective window
            choose=( object@dat@col.info$win==window | object@dat@col.info$win==window+1) # relates to the data matrix in object
            sum(choose)
            ttpos=object@dat@col.info$pos[choose]
            #
            # data matrix in object@dat may contain some SNPs in between that were not included in any large cluster
            choose2= (object@dat@col.info$pos >= min(ttpos) &  object@dat@col.info$pos <= max(ttpos))# relates to the data matrix in object@dat
            aapos=object@dat@col.info$pos[choose2]
            #aadat=object@dat@lib.freqs[choose2,]
            aadat=object@dat@lib.freqs[choose2,object@dat@use.libs,with=F]
            rownames(aadat)<-aapos
            
            dim(aadat)
            length(aapos)
            
            if (object@transf) {aadat=sqrt(aadat)}#; print("trans")}
            if (object@scaleSNP) {aadat=data.table(t(scale(t(aadat))))}#; print("scale")}
            
            a=nrow(aadat)
            aapos=aapos[rowSums(is.na(aadat))==0]
            aadat=aadat[rowSums(is.na(aadat))==0,] #na.omit(aadat)
            b=nrow(aadat)
            if (b<a) warning(paste("Due to NA positions (position of too low coverage)",round((1-b/a)*100,2),"% of the
                                   SNP positions were removed for the PCA analysis!"))
            
            # PCA
            pr.comp=prcomp(aadat,center = T, scale = T)
            
            # variances, i.e. proportion of variances explained
            vars <- pr.comp$sdev^2
            vars <- vars/sum(vars)
            
            # color according to performed clustering
            col=rep("black", length(aapos))
            # block ids for window
            hbs=rev_map(object)[names(rev_map(object))==paste0("w",window)][[1]]
            cc=rainbow(length(hbs))
            for (i in hbs)
            {
              col[ aapos %in% markers(object,i) ]=cc[which(hbs==i)]
            }
            
            # plot data points in terms of PCs
            par(mfrow=c(1,1), mar=c(4,4,1,1), oma=c(1,1,1,1))
            plot(pr.comp$x[,1], pr.comp$x[,2], pch=20, col=col
                 , xlab=paste("PC1,",round(vars[1],3)*100,"%")
                 , ylab=paste("PC2,",round(vars[2],3)*100,"%"))
            xx=par('usr') # figure x,ylim
            legend(0.85*(xx[2]-xx[1])+xx[1], 0.95*(xx[4]-xx[3])+xx[3]
                   , paste("hb",hbs), col=cc, pch=1, bty = "n")
          })


#' Hierachical clustering of haplotype-block marker for a window
#' 
#' Performs average linkage clustering for all SNPs in a window 
#' after filtering (i.e. \code{\link{initialize_SNP_time_series}}). 
#' 
#' @aliases inspect_window_avLink inspect_window_avLink.hbr
#' @author Susanne U. Franssen
#' @details Performs average linkage clustering on the time series data 
#' of the filtered SNPs and the selected time points (filtered and selected 
#' through \code{use.libs} in \code{\link{initialize_SNP_time_series}})
#' in a given window. The horizontal red line indicates the correlation
#' (1-correlation) threshold on which clusters were originally identified. 
#' Note: Clusters with less than \code{min.cl.cor} markers are discarded.
#' If no horizontal line is shown this indicates that all SNPs were assigned
#' to the same cluster.
#' @export
#' @importFrom stats hclust as.dist na.omit
#' @param object object of class \code{\link{hbr}} containing the results of 
#' reconstructed haplotype-blocks
#' @param window number of the window, which should be inspected
#' @param plotDendro Boolean indicating if dendrogram using average
#' linkage clustering should be plotted.
#' @param plotCluster Boolean indicating if clusters identified with average
#' linkage clustering should be visualized in PCs.
#' @param min.cl.cor minimum cluster correlation, if none is provided, the min.cl.cor
#' of the provided \code{\link{hbr}} object is taken.
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}}
#' \code{\link{inspect_window_dbScan}} \code{\link{number_hbr}}
#'
setMethod("inspect_window_avLink", "hbr", 
          function(object, window, min.cl.cor=0, plotDendro=T, plotCluster=T)
          {
            if (min.cl.cor==0){
              min.cl.cor=object@min.cl.cor
            }

            # extract data matrix from the respective window
            choose=( object@dat@col.info$win==window | object@dat@col.info$win==window+1) # relates to the data matrix in object
            sum(choose)
            ttpos=object@dat@col.info$pos[choose]
            #
            # data matrix in object@dat may contain some SNPs in between that were not included in any large cluster
            choose2= (object@dat@col.info$pos >= min(ttpos) &  object@dat@col.info$pos <= max(ttpos))# relates to the data matrix in object@dat
            aapos=object@dat@col.info$pos[choose2]
            #aadat=object@dat@lib.freqs[choose2,]
            aadat=object@dat@lib.freqs[choose2,object@dat@use.libs,with=F]
            rownames(aadat)<-aapos
            
            dim(aadat)
            length(aapos)
            
            if (object@transf) {aadat=sqrt(aadat)}#; print("trans")}
            if (object@scaleSNP) {aadat=data.table(t(scale(t(aadat))))}#; print("scale")}
            
            dsub.cor=cor(t(aadat),use="pairwise.complete.obs")
            dsub.cordist=1-as.dist(dsub.cor) # convert to dist format
            if (object@pos.cor) {dwin.cor[dwin.cor<0]=0}
            dsub.clust=hclust(dsub.cordist,method="average")
            dsub.clust$labels<-aapos
            
         
            # visulalize dendrogram
            if (plotDendro)
            {
              par(mfrow=c(1,1),mar=c(2,4,2.5,1))
              plot(dsub.clust, main="Average linkage clustering", ylab="Distance [1-correlation]", xlab="")
              abline(h=1-min.cl.cor, col="red")
            }
            if (plotCluster)
            {
              # visualize clusters on PCs
              dsub.cl=cutree(dsub.clust,h = 1-min.cl.cor)
              counts=table(dsub.cl) # count SNPs in each cluster
              small=as.numeric(names(counts))[counts<object@min.cl.size] # clusters that do not meet the minimum size
              dsub.cl[dsub.cl %in% small]<-0
              dsub.cl=dsub.cl+1
              #
              a=nrow(aadat)
              dsub.cl=dsub.cl[rowSums(is.na(aadat))==0]
              aadat=aadat[rowSums(is.na(aadat))==0,] #na.omit(aadat)
              b=nrow(aadat)
              if (b<a) warning(paste("Due to NA positions (position of too low coverage)",round((1-b/a)*100,2),"% of the
                                    SNP positions were removed for the PCA analysis! This can cause discrepancies in 
                                    the dendrogram and the PCA visualization!"))
              pr.comp=prcomp(aadat,center = T, scale = T, na.action=na.omit)
              size=rep(1.5,length(dsub.cl)) # point size: large are clusters small are discarded
              size[dsub.cl==1]=0.7
              comp=data.frame(pr.comp$x[,1:3]) 
              plot(comp, col=dsub.cl, pch=20, cex=size)
            }
            
          })



#' Clustering with the dbscan algorithm
#' 
#' Performs clustering with the dbscan (Density-based spatial clustering 
#' of applications with noise) for all SNPs in a window 
#' after filtering (i.e. \code{\link{initialize_SNP_time_series}})
#' and visualizes found clusters based on principal components. 
#' NOTE: The clusters identified here are not nessesarily identical
#' with the clusters identified with average linkage clustering. Therefore
#' if haplotype reconstruction was done using average linkage clustering
#' the clusters shown here can be different from clusters of identified haplotype 
#' blocks. 
#' 
#' @aliases inspect_window_dbScan inspect_window_dbScan.hbr
#' @author Susanne U. Franssen
#' @details Performs clustering with \href{https://CRAN.R-project.org/package=dbscan}{dbscan} 
#' (Density-based spatial clustering of applications with noise) for all SNPs in a window 
#' after filtering (i.e. \code{\link{initialize_SNP_time_series}})
#' and visualizes found clusters based on principal components. 
#' NOTE: The clusters identified here are not nessesarily identical
#' with the clusters identified with average linkage clustering. Therefore
#' if haplotype reconstruction was done using average linkage clustering
#' the clusters shown here can be different from clusters of identified haplotype 
#' blocks. This method is rather indendet for inspection and getting an idea
#' if haplotype reconstruction should rather be run using dbscan instead 
#' of average linkage clustering.
#' 
#' This package used the dbscan implementation of the package \href{https://CRAN.R-project.org/package=dbscan}{dbscan}
#' (Michael Hahsler) originally described by Ester et al. (1996).
#' @export
#' @importFrom dbscan dbscan kNNdistplot
#' @param object object of class \code{\link{hbr}} containing the results of 
#' reconstructed haplotype-blocks
#' @param window Number of the window, which should be inspected
#' @param minPts Number of minimum points in the eps region (for core points).
#' (default \code{min.cl.size} used for haplotype block reconstruction)
#' @param eps The size of the epsilon neighborhood. Defines the distance between
#' samples that built up a core cluster, for details see ?\code{dbscan} of
#' the \code{dbscan} package. (This value is only used for the \code{plotCluster} 
#' functionality.)
#' @param plotkNN Boolean indicating if the kNN distance plot should be printed
#' , for details see ?\code{kNNdistplot} of the \code{dbscan} package. Briefly,
#' the y-value where a bent in the curve is visible is a good indicator of the eps
#' value to choose for the given k = minPts.
#' @param plotCluster Boolean indicating if clusters identified with average
#' linkage clustering should be visualized in PCs.
## @references Martin Ester, Hans-Peter Kriegel, Joerg Sander, Xiaowei Xu (1996). 
## A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases 
## with Noise. Institute for Computer Science, University of Munich. Proceedings 
## of 2nd International Conference on Knowledge Discovery and Data Mining (KDD-96).
#' @seealso \code{\link{hbr}} \code{\link{ex_dat}} \code{\link{summary.hbr}} \code{\link{plot.hbr}} 
#' \code{\link{plot_cluster_trajectories}} \code{\link{plot_marker_trajectories}}
#' \code{\link{map}} \code{\link{rev_map}} \code{\link{markers}} \code{\link{plot_hbr_freq}}
#' \code{\link{inspect_window}} \code{\link{inspect_window_PCA}} \code{\link{inspect_window_avLink}}
#' \code{\link{number_hbr}}
#'
setMethod("inspect_window_dbScan", "hbr", 
          function(object, window, minPts=object@min.cl.size, eps, plotkNN=T, plotCluster=T)
          {
            # extract data matrix from the respective window
            choose=( object@dat@col.info$win==window | object@dat@col.info$win==window+1) # relates to the data matrix in object
            sum(choose)
            ttpos=object@dat@col.info$pos[choose]
            #
            # data matrix in object@dat may contain some SNPs in between that were not included in any large cluster
            choose2= (object@dat@col.info$pos >= min(ttpos) &  object@dat@col.info$pos <= max(ttpos))# relates to the data matrix in object@dat
            aapos=object@dat@col.info$pos[choose2]
            #aadat=object@dat@lib.freqs[choose2,]
            aadat=object@dat@lib.freqs[choose2,object@dat@use.libs,with=F]
            rownames(aadat)<-aapos
            
            dim(aadat)
            length(aapos)
            
            if (object@transf) {aadat=sqrt(aadat)}#; print("trans")}
            if (object@scaleSNP) {aadat=data.table(t(scale(t(aadat))))}#; print("scale")}
            
            a=nrow(aadat)
            aapos=aapos[rowSums(is.na(aadat))==0]
            aadat=aadat[rowSums(is.na(aadat))==0,] #na.omit(aadat)
            b=nrow(aadat)
            if (b<a) warning(paste("Due to NA positions (position of too low coverage)",round((1-b/a)*100,2),"% of the
                                   SNP positions were removed for the PCA analysis!"))
            
            if (plotkNN){
              kNNdistplot(aadat, k=minPts)
            }
            
            if (plotCluster){
              xx=dbscan(aadat, eps=eps, minPts = minPts)
              #
              pr.comp=prcomp(aadat,center = T, scale = T)
              size=rep(1.5,length(xx$cluster)) # point size: large are clusters small are discarded
              size[xx$cluster==0]=0.7
              comp=data.frame(pr.comp$x[,1:3]) 
              plot(comp, col=xx$cluster+1, pch=20, cex=size)
            }
          })




