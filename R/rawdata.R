#' Import .mzXML data into a list
#'
#' \code{rawdata} reads a single .mzXML file into a list object which consists of a named list of data.table objects. These lists contain the raw mass spectral data from the file.
#'
#' Also handles scan stitching if multiple scan segments were employed.
#'
#' @param file Character. Path to the .mzXML file
#' @param splits Matrix with N rows and 2 columns containing the mass ranges to extract from sequential scans.  Column 1 is the start of the mass range and column 2 is the end of the mass range.  The first row will apply to the first scan in the data file, the second to the second, and so on.  When the matrix rows are exhausted it will return to the first row. For example, a two row matrix, the first row will apply to all odd scans and the second to all even scans. Each set of N scans will be combined into a single scan for further analysis.
#'
#' @return List containing lists named k (mass peaks), s (scan information), and metadata.
#'
#' @examples
#' \dontrun{
#' rawdata('file1.mzXML', splits = rbind(c(0, 300), c(300, 3000)))
#' }
#'
#' @export
#'
rawdata = function(file, splits = NULL) {
  cat("\n\nLoading and merging raw data from", file, "\n")

  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Package 'mzR' needed for this function to work. Please install it.", call. = FALSE)
  }

  suppressWarnings(suppressMessages(library(mzR)))

  ramp = mzR:::rampOpen(file)
  data = mzR:::rampRawData(ramp)
  mzR:::rampClose(ramp)

  if (!is.null(splits)) {
    nsplits = nrow(splits)

    if (nsplits > 1) { message("Scan stitching was applied! Scan numbers no longer corrspond to original file.") }

    seglengths = diff(c(data$scanindex, length(data$mz)))
    seg = rep(rep(seq(nsplits), ceiling(length(data$rt)/nsplits))[1:length(seglengths)], seglengths)
    keepmat = outer(splits[,1], data$mz, "<") & outer(splits[,2], data$mz, ">") & outer(seq(nsplits), seg, "==")
    keepvec = matrixStats::colAnys(keepmat)

    newindices = cumsum(keepvec)[data$scanindex+1]-1
    newindices = newindices[seq(from = 1, to = length(data$scanindex), by = nsplits)]

    newis = seq(from = 1, to = length(data$rt), by = nsplits)

    data$mz = data$mz[keepvec]
    data$intensity = data$intensity[keepvec]
    data$scanindex = as.integer(newindices)
    data$rt = data$rt[newis]
    data$acquisitionNum = seq_along(data$acquisitionNum[newis])
    data$polarity = as.integer(data$polarity[newis])
  }

  data$acquisitionNum = seq_along(data$acquisitionNum)

  list(
    s = data.table(data.frame(s = data$acquisitionNum, rt = data$rt, polarity = data$polarity)),
    k = data.table(data.frame(mz = data$mz, i = data$intensity, s = rep(data$acquisitionNum, diff(c(data$scanindex, length(data$mz)))), k = seq_along(data$mz))),
    metadata = list(file=file)
  )
}
