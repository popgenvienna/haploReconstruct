
# checks validity of the input signature for plot_hbr_freq method
#
validity_plot_hbr_freq<-function(object, hbr_id, replicate, timepoint, window
                                 , cols, add, cex
                                 , xlab, ylab, xlim, ylim)
{
  # parameter presence
  if (is.null(object)|is.null(hbr_id)|is.null(replicate)|is.null(timepoint)|is.null(window)
      |is.null(xlab)|is.null(ylab)|is.null(ylim)|is.null(cex))
  { stop("Parameter missing with no default.") }
  
  # hbr_id
  if (!is.numeric(hbr_id)){
    stop(paste("Please provide an integer in the range from 1 to",
               length(object@cl.long.m),"for hbr_id."))
  } else {
    if(!hbr_id%%1==0) stop(paste("Please provide an integer in the range from 1 to",
                                 length(object@cl.long.m),"for hbr_id."))
  }
  
  if (!is.logical(add)) stop("add needs to be of type logical.")
  
}



# used in show
#
print_wincl_hbr_info=function(ll)
{
  v=unlist(lapply(ll,length))
  m=matrix(c(1:length(v),v),ncol=2)
  m=m[m[,2]!=0,]
  for (r in 1:nrow(m))
  {
    cat("\t",m[r,1],"\t",m[r,2],"\n")
  }
}



