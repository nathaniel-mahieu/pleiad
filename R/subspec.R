#' Return nearest mass peaks to supplied peaks in a set of spectra.
#'
#' \code{specs_nearestmasses} for every spectra returns the peak nearest the supplied mass.
#'
#' @param masses data.table containing column mz
#' @param specs data.table containing columns mz, s
#'
#' @return data.table containing columns mz.target, mz, ppm, s, as well as original columns.
#'
#' @export
#'
subspec_nearestmasses = function(masses, specs) {
  masses[,mz.target := mz]
  specs[,mzs := mz]
  setkey(specs, mz)

  tmp = specs[,.SD[masses, ,roll="nearest", on=.(mz)], by=.(s)][,mz:=mzs][,!"mzs"]
  tmp[,ppm := (mz - mz.target)/mz * 1E6]
  tmp[,subspec:=masses$mz %>% round(4) %>% paste(collapse="x") %>% factor]
  tmp[,i.norm := i/max(i),by=.(s)]

  tmp
}

#' Correct deviation from target masses within subspectra
#'
#' \code{subspec_masscorrect} for every subspec corrects the average mz.target - mz
#'
#' @param specs data.table containing columns mz, subspec
#' @param maxppm the maximum ppm error for an observation to be used in calculating the correction
#'
#' @return data.table containing new column mz.observed and a corrected mz value
#'
#' @export
#'
subspec_masscorrect = function(specs, maxppm = Inf) {
  specs[,mz.observed := mz]
  specs[,ppm.observed := ppm]

  corrections = specs[align.rtuserange == T][abs(ppm) < maxppm][,.(mzd = mean(mz.target-mz)),by=.(subspec)]

  specs[corrections,mz := mz + mzd,on=.(subspec)]
  specs[,ppm := (mz - mz.target) / mz.target * 1E6]

  specs
  }


#' Find a the average isotopologue spectrum considering all spectra in the specified range.
#'
#' \code{subspec_spectra} for every subspec finds the mean and standard deviation of each observed isotopologue.
#'
#' @param specs data.table containing columns mz, subspec, align.rtuserange, i, s, ppm, rt, mz.target
#' @param maxppm the maximum ppm error for an observation to be used in calculating the correction
#' @param s.aroundmax the number of scans on either side of the maxima in the specified range to aggregate. The default, Inf, takes the whole specified range.
#' @param id.cols columns to carry through
#'
#' @return data.table containing each averaged spectra, including normalized intensity on a scan to scan basis
#'
#' @export
#'
subspec_spectra = function(specs, maxppm = Inf, s.aroundmax = Inf, id.cols = NULL) {

  s.peak = specs[align.rtuserange == T][,.(i.tot = sum(i)),by=.(s)][which.max(i.tot)]$s
  specs[,i.norm := i/max(i),by=.(s)]

  tmp = specs[abs(ppm) < maxppm][align.rtuserange == T][abs(s-s.peak) < s.aroundmax][
    ,.(
      mz = mean(mz),
      rt = mean(rt),
      ppm = mean(ppm),
      i = mean(i),
      i.sd = sd(i),
      i.norm=mean(i.norm, na.rm=T),
      i.norm.sd = sd(i.norm, na.rm=T)
    ),by=.(subspec, mz.target)]

  tmp.extracols = specs[abs(ppm) < maxppm][align.rtuserange == T][abs(s-s.peak) < s.aroundmax][,.SD[,id.cols,with=F][1],by=.(subspec, mz.target)]

  tmp[tmp.extracols,,on=.(subspec, mz.target)]
}
