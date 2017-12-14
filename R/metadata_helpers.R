make_metadata_common_name = function(files) {

  x.l = substr(basename(files), start = 1, stop = nchar(basename(files))-6) %>% strsplit(split="_") # Split by "_"

  mat = sapply(x.l, function(x) { x == x.l[[1]] }) %>% do.call (what=rbind)

  keeps = which(colSums(mat) == nrow(mat))

  if (keeps < 1) { stop("No name could be extracted.") }

  paste(x.l[[1]][keeps], collapse="_")
}

make_metadata_csv_skeleton = function(files) {
  x.l = substr(basename(files), start = 1, stop = nchar(basename(files))-6) %>% strsplit(split="_") # Split by "_"

  n = sapply(x.l, length) %>% max
  x.l = lapply(x.l, function(x) { length(x) = n; x })

  x.df = x.l  %>% do.call(what=rbind)

  colnames(x.df) = paste(sep=".", "NAME.ME", seq_len(ncol(x.df)))
  x.df = cbind(x.df, file = files)

  x.df %>% write.csv(file="metadata_csv_skeleton.csv", row.names = F)

  x.df
}