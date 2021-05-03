
library("flowCore")
library("flowViz")
library("ggcyto")

frames <- lapply(dir("data/flow_cytometry/FCS/", full.names = TRUE), read.FCS)
names(frames) <- sapply(frames, keyword, "TUBE NAME")
fs <- as(frames, "flowSet")
sampleNames(fs)

phenoData(fs)

fsApply(fs, each_col, sum)

f <- fs[[10]]
f
autoplot(transform(f, 
                   "Pacific Blue-H" = log1p(`Pacific Blue-H`), 
                   "APC-H" = log1p(`APC-H`)),
         "Pacific Blue-H", "APC-H")

?autoplot         
