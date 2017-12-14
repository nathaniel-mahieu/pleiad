lassign = function(..., pos = 1, envir = as.environment(pos)) {
  as = list(...)

  for (i in seq_along(as)) {
    assign(names(as[i]), as[[i]], envir = envir)
  }
}
