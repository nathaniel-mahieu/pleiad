# Update if necessary
# setRepositories(ind = 1:2); devtools::install_github("nathaniel-mahieu/pleiad")

# Setup
setwd()
library(pleiad) # Load the package
source("workflow_functions.R")
theme_set(theme_light() + theme(strip.background = element_rect(fill="white"), strip.text = element_text(colour="black"), legend.position="top")) # Make ggplot plots look pretty



# Configure parallel cluster
library(doParallel)
cl = makeCluster(3) # 3 is the number of cores
registerDoParallel(cl)


# Load and test metabolite library
fread("compounds_nate.csv") %>% assign(x="metabolite_dt", envir=pleiad.env)

# Examples
m("glu")
mz("glu", .ion = "+Na")
mass("glu")
isotopologues("glu", .atom = "C", .ion = "-H")[1:6]
isotopologues.combn(mz("glu"), .ion = "-H", C=0:4, N=0:2, H=0) %>% print


# Files to process
# First step in the analysis is to load all of the data files into R
files = list.files("PATH/TO/FILES", pattern = "*\\.mzXML$", full.names = T)
files

make_metadata_csv_skeleton(files) # Go find, copy, and edit "metadata_csv_skeleton.csv" this making sure the metadata is named and accurate


filename_to_metadata_vector = function(file.) { # you can replace this code with anything that takes a filename and returns the metadata you want associated with the sample as a named list.
  metadata.dt = read.csv("metadata_csv_actual.csv") %>% data.table # Make this point to your metadata file.

  metadata.dt[file == file.] %>% c
  }


ms_list.l = foreach(file = files, .packages = "pleiad") %dopar% {
  data = rawdata(file, splits = c(0,400,400,3000) %>% matrix(ncol = 2, byrow = T)) # ADJUST: splits =

  data$metadata$factors = filename_to_metadata_vector(file)

  data
  }


###########
# Analyze data with the metabolite helpers
timestamp. = round(as.numeric(Sys.time()))
metab = m("glu")
metab

# Pick the file on which the suppled retention times are based
align.file = ms_list.l[[1]]

# Grab the mz of interest
mz. = mz("glu", .ion = "-H")
align.trace = align.file$k[,.SD[which.min(abs(mz-mz.))],by=.(s)][align.file$s,,on=.(s)][,.(mz, i, rt)]
align.trace[abs(rt-200)<100,.(rt, i)][,.(rt,i)] %>%  ggplot() + geom_line(aes(x=rt, y = i))

# Extract the data
analysis = targeted_isotopologues(ms_list.l = ms_list.l, target.mz = mz., align.trace = align.trace, align.window.size=10, align.rtrange = c(200, 250), align.rtuserange = c(220, 240), align.saroundmax = 6, C =0:6, N = 0:2, H = 0, maxppm = 4)

summarize_targeted_isotopologues(analysis, name = metab$name, timestamp = timestamp.)
plot_targeted_isotopologues(analysis$spectra)
plot_targeted_isotopologues_qc(analysis)
