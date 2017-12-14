#' Convenience function for data processing
#'
#' \code{dtw_getspecs} aligns a mass trace to a supplied trace and returns a specified retention time range post alignment.
#'
#' @param ms_list.l A list of ms_lists as obtained by \code{rawdata}
#' @param metadata A data.table of metadata with column file corresponding to the index of ms_list.l
#' @param target.mz The M0 peak mass
#' @param align.trace A data.table of columns rt and i of the EIC which to align all files to
#' @param align.rtrange Numeric vector containing the retention time range which to align. c(rtmin, rtmax)
#' @param align.rtuserange Numeric vector containing the retention time range which, after aligned, to use in calculating abundances. c(rtmin, rtmax)
#' @param align.saroundmax Numeric of the number of scans before and after the maximum in the suppled align.rtuserange which to take for analysis. If Inf (default) takes the whole supplied region.  Notably, this could cause the script to grab a nearby, large peak if it is in the region. This seems rare but we can implement better logic if necessary.
#' @param C.H.N Numeric vectors containing the number of isotopologues for each atom to consider. Ex. 0:3 considers M0, M1, M2, and M3.
#' @param maxppm Numeric the maximum ppm within which to consider a peak
#'
#' @return List containing spectra (averaged/combined) and subspecs (original peaks per scans in align.rtrange)
#'
#' @description 1. Aligns a mass trace from the file to be analyzed with a supplied trace such that the same region can be identified between the two. 2. From the desired region the closest peak to each suppled mass is grabbed - desired masses are defined by C, N, H and target.mz. 3. Spectra are normalized to the M0 and mass corrected.  4. Spectra are averaged to provide an isotopologue spectrum for the species.  Note: if mass error is too high the wrong peaks will be grabbed.  Incorrect peaks within maxppm will be returned and plotted as real peaks. The rt region selection logic is robust but should be checked for errors for each analyte.
#'
#' @export
#'
targeted_isotopologues = function(ms_list.l, target.mz, align.trace, align.rtrange, align.rtuserange, align.saroundmax = Inf, align.window.size, C = 0, N = 0, H = 0, maxppm = 4) {

  anal.l = foreach(ms_list = ms_list.l, .packages = c("pleiad")) %dopar% { cat(ms_list$metadata$file, " | ")
    anal = dtw_getspecs(ms_list = ms_list, align.mz = target.mz, align.trace = align.trace, align.rtrange = align.rtrange, align.rtuserange = align.rtuserange, window.size = align.window.size)

    anal$subspecs = subspec_nearestmasses(isotopologues.combn(target.mz, C=C, N=N, H=H), anal$specs)
    anal$subspecs = subspec_masscorrect(anal$subspecs, maxppm = maxppm)

    anal$spectra = subspec_spectra(anal$subspecs, s.aroundmax = align.saroundmax, id.cols=c("C", "N", "H", "align.rtuserange", "align.rtrange"), maxppm = maxppm)

    anal$metadata = ms_list$metadata$factors

    anal
  }

  tmp = list()
  tmp$rawdata = anal.l

  tmp$spectra = lapply(seq_along(anal.l), function(i.) {
    anal.l[[i.]]$spectra[,file:=i.]
  }) %>% do.call(what=rbind)

  metadata = lapply(anal.l, '[[', "metadata") %>% do.call(what = rbind) %>% data.table %>% { colnames(.) = names(anal.l[[1]]$metadata); .[,file := seq_len(.N)] }

  tmp$spectra = tmp$spectra[metadata,,on=.(file),nomatch=0]

  tmp
}


plot_targeted_isotopologues_qc = function(analysis) {
  # Region QC
  traces = lapply(seq_along(analysis$rawdata), function(i.) {
    x = analysis$rawdata[[i.]]

    rtrange = x$specs$rt %>% range(na.rm = T)
    x$target.trace[,inexpandedrange := findInterval(rt, rtrange + c(-10, 10)) == 1]

    x$target.trace[,':='(file=i.)]
  }) %>% do.call(what = rbind)

  { traces[inexpandedrange==T] %>% ggplot() + geom_line(aes(x = rt, y = i, colour = interaction(align.rtrange, align.rtuserange), group = file)) +
    facet_wrap(~file, scales = "free_y") + #theme(legend.position="none") +
    theme(axis.text.x=element_blank(), axis.text.y=element_blank()) +
    ggtitle("QC: Regions utilized for analysis", "Peak region detection may make mistakes - double check the appropriate region was considered if you see an interesting difference.") } %>% print


  Sys.sleep(2)

  # Trap fill QC
  traces = lapply(seq_along(analysis$rawdata), function(i.) {
    x = analysis$rawdata[[i.]]
    y = analysis$rawdata[[i.]]$specs
    if (nrow(y) < 1) {return(NULL)}

    rtrange = x$specs$rt %>% range(na.rm = T)

    z = y[,quantile(i, c(0.0025, .25)) %>% as.list,by=.(rt)]
    z[,':='(file=i.)]
    z[,inexpandedrange := findInterval(rt, rtrange + c(-10, 10)) == 1]

  }) %>% do.call(what = rbind) %>% melt(id.vars=c("rt", "file", "inexpandedrange"))

  { traces[inexpandedrange==T] %>% ggplot() + geom_line(aes(x = rt, y = value, colour = variable, group = variable)) +
    facet_wrap(~file, scales = "free_y") + theme(legend.position="top") +
    #theme(axis.text.x=element_blank(), axis.text.y=element_blank()) +
    ggtitle("QC: Peak Intensity Quantiles. Approximates limit of detection.", "Indicator of trap filling. Some small abundance peaks may be harmonics.") } %>% print

  Sys.sleep(2)

  # Mass QC
  analysis$spectra[,mz.target := round(mz.target, 5)]
  foo = analysis$spectra[,.(mean = mean(mz), ppm5 = mean(mz)/1E6*5, meanrt = mean(rt), rt2 = 2),by=.(mz.target)]

  { analysis$spectra %>% ggplot() + geom_point(aes(x = mz, y = rt, colour = factor(file))) + facet_wrap(~mz.target, scales = "free") + ggtitle("QC: Masses used for statistics", "Cross indicates a 2.5/5 ppm in x and 2 seconds in y. Large drift might indicate the incorrect peak was used.") +
    geom_segment(data = foo, aes(x=mean - ppm5/2, xend = mean + ppm5/2, y = meanrt, yend = meanrt), colour = "grey")+
    geom_segment(data = foo, aes(x=mean - ppm5/2/2, xend = mean + ppm5/2/2, y = meanrt, yend = meanrt), colour = "darkgrey")+
    geom_segment(data = foo, aes(x=mean, xend = mean, y = meanrt-rt2/2, yend = meanrt+rt2/2), colour = "grey") } %>% print


}
