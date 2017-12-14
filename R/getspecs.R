#' Get spectra after DTW alignment
#'
#' \code{dtw_getspecs} aligns a mass trace to a supplied trace and returns a specified retention time range post alignment.
#'
#' @param ms_list Raw data list containing list items k and s.
#' @param align.mz The mass in ms_list to align to.
#' @param align.trace Trace data.table containing columns rt and i to align to.
#' @param align.rtrange Vector of length 2 with begenning and end of retention time range to return. Range referse to the trace in align.trace.
#'
#' @return List containing specs, rt_map, and params.
#'
#' @export
#'
dtw_getspecs = function(ms_list, align.mz, align.trace, align.rtrange, align.rtuserange, window.size) {
  target.trace = ms_list$k[,.SD[which.min(abs(mz-align.mz))],by=.(s)][ms_list$s,,on=.(s)][,.(mz, i, rt)]

  setorder(target.trace, rt)
  setorder(align.trace, rt)

  rt_map = {

    rts = range(c(align.trace$rt, target.trace$rt)) %>% { seq(.[1], .[2], mean(diff(align.trace$rt))/2) }

    rts = rts[findInterval(rts, align.rtrange %>% { . + c(-1, 1)*window.size*2 }) == 1]

    target.trace.inter = approx(target.trace[,.(rt, i)], xout=rts)$y %>% filter(1/rep(5,5)) %>% { .[is.na(.)] = 0; . } %>% { ./max(., na.rm=T) }
    align.trace.inter = approx(align.trace[,.(rt, i)], xout=rts)$y %>% filter(1/rep(5,5)) %>% { .[is.na(.)] = 0; . } %>% { ./max(., na.rm=T) }

    dt = dtw(target.trace.inter, align.trace.inter,
             open.end = F, open.begin = F, window = "sakoechiba", window.size = ceiling(window.size/mean(diff(rts))),
             step.pattern = symmetric2, keep.internals = F
    )

    data.table(target.rt = rts[dt$index1], align.rt = rts[dt$index2])
  }

  target.rtrange = rt_map[findInterval(align.rt, align.rtrange) == 1]$target.rt %>% range
  target.align.rtuserange = rt_map[findInterval(align.rt, align.rtuserange) == 1]$target.rt %>% range

  target.specs = ms_list$k[ms_list$s,,on=.(s)][findInterval(rt, target.rtrange)==1]
  target.specs[,':='(align.rtrange = T, align.rtuserange = findInterval(rt, target.align.rtuserange)==1)]
  target.trace[,':='(align.rtrange = findInterval(rt, target.align.rtuserange)==1, align.rtuserange = findInterval(rt, target.rtrange)==1)]


  list(
    specs = target.specs,
    rt_map = rt_map,
    align.trace = align.trace,
    target.trace = target.trace,
    params = list(align.mz = align.mz, align.trace = align.trace, align.rtrange = align.rtrange)
  )

}