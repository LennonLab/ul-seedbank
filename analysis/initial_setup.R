rm(list = ls())
#source("bin/mothur_tools.R")
require(vegan)
require(reshape2)
require(tidyverse)
require(lubridate)

# OTUs <- read.otu("data/ul.bac.final.shared")
# OTUs.tax <- read.tax("data/ul.bac.final.0.03.taxonomy")

# Make the design matrix
# design <- data.frame(sample.name = rownames(OTUs))
# design <- cbind.data.frame(design,
#                 colsplit(design$sample.name, pattern = "[c D]", c("sample.type", "sample.id")))
# design$sample.type <- substr(design$sample.name, start = 3, 3)
# for(each in 1:length(design$sample.type)){
#   if(design$sample.type[each] == "D") design$sample.type[each] <- "DNA"
#   if(design$sample.type[each] == "c") design$sample.type[each] <- "RNA"
# }

# Sort by sample type and number
# OTUs.sort <- cbind(design, OTUs)
# OTUs.sort <- arrange(OTUs.sort, sample.type, sample.number)
# design <- arrange(design, sample.type, sample.number)
# OTUs.sort <- column_to_rownames(OTUs.sort, var = "sample.ID")
# OTUs <- OTUs.sort[,-c(1:2)]
# 
# write.csv(design, file = "data/design.csv", row.names = F)
# saveRDS(OTUs, file = "data/Rdata/OTUs.rda")
# saveRDS(OTUs.tax, file = "data/Rdata/OTUtax.rda")
OTUs.tax <- readRDS(file = "data/Rdata/OTUtax.rda")
OTUs <- readRDS(file = "data/Rdata/OTUs.rda")
design <- read.csv(file = "data/design.csv")

# Make rel abund matrices and split into active total comms
OTUs <- OTUs[,-which(colSums(OTUs) < 3)]

OTUs <- rrarefy(OTUs, sample = min(rowSums(OTUs)))

OTUs.active <- OTUs[which(design$sample.type == "RNA"),]
OTUs.total <- OTUs[which(design$sample.type == "DNA"),]
#OTUs[which(design$sample.type == "DNA"),] <- decostand(OTUs.active, method = "pa") + OTUs.total

OTUs.REL <- decostand(OTUs, method = "hellinger")
OTUs.REL.total <- OTUs.REL[which(design$sample.type == "DNA"),]
OTUs.REL.active <- OTUs.REL[which(design$sample.type == "RNA"),]

# read in environmental data
env.data <- read.table("data/ul-seedbank.env.txt", sep="\t", header=TRUE)
env.data$date <- parse_date_time(env.data$date, "m d y")
