pleiad.env = new.env(parent = emptyenv())
assign("metabolite_dt", data.table(name = c("glutamate", "aspartate"), abbreviation = c("glu","asp"), formula = c("C5H9NO4", "C4H7NO4")), envir=pleiad.env)
assign("metabolite_ions", data.table(ion = c("-H", "+H", "+Na"), shift = c(-1.007825, +1.007825, 22.98977)), envir=pleiad.env)
assign("metabolite_shifts", data.table(atom = c("C", "H", "N"), shift = c(1.003355, 1.006277, 0.9970349)), envir=pleiad.env)

#' Retrieve an entry from the metabolite_dt
#'
#' \code{m} Searches metabolite_dt in environment pleiad.env for rows which match supplied search either by name or abbreviation.
#'
#' @param search Character. Name or abbreviation string.  Case sensitive.
#'
#' @return data.table containing a single row from metabolite_dt
#'
#' @export
#'
m = function(search) {
  matches = get("metabolite_dt", envir=pleiad.env)[,c(which(name == search), which(abbreviation == search))]
  if (length(matches) > 1) warning("Search matched multiple metabolites.")
  if (length(matches) < 1) warning("Search matched no metabolites.")

  get("metabolite_dt", envir=pleiad.env)[matches[1]]
}

#' Retrieve a single column from the metabolite_dt
#'
#' \code{m} Searches metabolite_dt in environment pleiad.env for rows which match supplied search either by name or abbreviation and returns the specified col.
#'
#' @param search Character. Name or abbreviation string.  Case sensitive.
#' @param col Character. Name of the desired column.
#'
#' @return value from metabolite_dt
#'
#' @export
#'
retrieve.metabolite_dt = function(search, col="formula") {
  m(search)[,col,with=F][[1]]
}

#' Retrieve a m/z from the metabolite_dt
#'
#' \code{m} Searches metabolite_dt in environment pleiad.env for rows which match supplied search either by name or abbreviation and returns the mass corresponding to that formula.
#'
#' @param search Character. Name or abbreviation string.  Case sensitive.
#' @param .ion Character. Name of the desired ion transformation. Must be defined in metabolite_ions
#'
#' @return m/z corresponding to the supplied formula and ion.
#'
#' @export
#'
mz = function(search, .ion = "-H") {
  shift = get("metabolite_ions", envir=pleiad.env)[ion == .ion, shift]

  if (length(shift) == 0) warning("Specified ion not found: ", .ion)

  mass(search) + shift
}

#' Retrieve a neutral mass from the metabolite_dt
#'
#' \code{m} Searches metabolite_dt in environment pleiad.env for rows which match supplied search either by name or abbreviation and returns the mass corresponding to that formula.
#'
#' @param search Character. Name or abbreviation string.  Case sensitive.
#'
#' @return mass corresponding to the supplied formula
#'
#' @export
#'
mass = function(search) {
 Rdisop::getMolecule(retrieve.metabolite_dt(search, "formula"))$exactmass
}

#' Calculate isotopologues given a metabolite from metabolite_dt, an ion, and the atom for which isotopes have been substituted.
#'
#' \code{m} Searches metabolite_dt in environment pleiad.env for rows which match supplied search either by name or abbreviation and returns the mass corresponding to that formula.
#'
#' @param search Character. Name or abbreviation string.  Case sensitive.
#' @param .ion Character. Name of the desired ion transformation. Must be defined in metabolite_ions
#' @param .atom Character. Name of the desired atom's isotope masses. Must be defined in metabolite_shifts
#' @param maxn Numeric. The largest isotopologue desired
#'
#' @return m/z vector corresponding to the supplied formula and ion and shift.
#'
#' @export
#'
isotopologues = function(search, .atom = "C", .ion="-H", maxn = 100) {
  shift = get("metabolite_shifts", envir=pleiad.env)[atom == .atom, shift]

  if (length(shift) == 0) warning("Specified atom not found: ", .atom)

  mz(search, .ion) + shift * 0:maxn
}

#' Calculate isotopologues given a metabolite from metabolite_dt, an ion, and the atom for which isotopes have been substituted.
#'
#' \code{m} Searches metabolite_dt in environment pleiad.env for rows which match supplied search either by name or abbreviation and returns the mass corresponding to that formula.
#'
#' @param search Character. Name or abbreviation string.  Case sensitive.
#' @param .ion Character. Name of the desired ion transformation. Must be defined in metabolite_ions
#' @param .atom Character. Name of the desired atom's isotope masses. Must be defined in metabolite_shifts
#' @param maxn Numeric. The largest isotopologue desired
#'
#' @return m/z vector corresponding to the supplied formula and ion and shift.
#'
#' @export
#'
isotopologues.combn = function(mz., .ion="-H", C=0:6, N=0:2, H=0) {
  metabolite_shifts = get("metabolite_shifts", envir=pleiad.env)

  atom_mat = expand.grid(C=C, N=N, H=H)
  shift_mat = matrix(rep(metabolite_shifts[match(atom, names(atom_mat)), shift], each=nrow(atom_mat)), ncol=ncol(atom_mat))

  shifts = rowSums(atom_mat * shift_mat)

  shift.dt = data.table(atom_mat)
  shift.dt[,shift:=shifts]

  shift.dt[ , mz := mz. + shift]

  shift.dt
}

count.elements <- function(formula) {
  # count the elements in a chemical formula   20120111 jmd
  # this function expects a simple formula,
  # no charge or parenthetical or suffixed subformulas
  # regular expressions inspired by an answer on
  # http://stackoverflow.com/questions/4116786/parsing-a-chemical-formula-from-a-string-in-c
  #elementRegex <- "([A-Z][a-z]*)([0-9]*)"
  elementSymbol <- "([A-Z][a-z]*)"
  # here, element coefficients can be signed (+ or -) and have a decimal point
  elementCoeff <- "((\\+|-|\\.|[0-9])*)"
  elementRegex <- paste(elementSymbol, elementCoeff, sep="")
  # stop if it doesn't look like a chemical formula
  validateRegex <- paste("^(", elementRegex, ")+$", sep="")
  if(length(grep(validateRegex, formula)) == 0)
    stop(paste("'",formula,"' is not a simple chemical formula", sep="", collapse="\n"))
  # where to put the output
  element <- character()
  count <- numeric()
  # from now use "f" for formula to make writing the code easier
  f <- formula
  # we want to find the starting positions of all the elemental symbols
  # make substrings, starting at every position in the formula
  fsub <- sapply(1:nchar(f), function(i) substr(f, i, nchar(f)))
  # get the numbers (positions) that start with an elemental symbol
  # i.e. an uppercase letter
  ielem <- grep("^[A-Z]", fsub)
  # for each elemental symbol, j is the position before the start of the next
  # symbol (or the position of the last character of the formula)
  jelem <- c(tail(ielem - 1, -1), nchar(f))
  # assemble the stuff: each symbol-coefficient combination
  ec <- sapply(seq_along(ielem), function(i) substr(f, ielem[i], jelem[i]))
  # get the individual element symbols and coefficients
  myelement <- gsub(elementCoeff, "", ec)
  mycount <- as.numeric(gsub(elementSymbol, "", ec))
  # any missing coefficients are unity
  mycount[is.na(mycount)] <- 1
  # append to the output
  element <- c(element, myelement)
  count <- c(count, mycount)
  # in case there are repeated elements, sum all of their counts
  # (tapply hint from https://stat.ethz.ch/pipermail/r-help/2011-January/265341.html)
  # use simplify=FALSE followed by unlist to get a vector, not array 20171005
  out <- unlist(tapply(count, element, sum, simplify=FALSE))
  # tapply returns alphabetical sorted list. keep the order appearing in the formula
  out <- out[match(unique(element), names(out))]

  out.zeros = c(C=0, N=0, H=0, P=0, O=0, S=0)

  return(c(out, out.zeros[!names(out.zeros) %in% names(out)]))
}
