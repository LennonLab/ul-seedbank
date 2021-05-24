
library("flowCore")
library("flowViz")
library("flowStats")
library("ggcyto")
library("tidyverse")
library("patchwork")
library("lubridate")
library("scales")

design <- read.csv(file = "data/design.csv")

env.data <- read.table("data/ul-seedbank.env.txt", sep="\t", header=TRUE)
env.data$date <- as.Date(parse_date_time(env.data$date, "m d y"))
env.data$doc[which(env.data$doc == "**")] <- NA
env.data$doc <- as.numeric(as.character(env.data$doc))
summary(env.data)
design <- left_join(design, env.data[,c("sample.id", "date")])
env.data <- env.data[-which(!(env.data$date %in% design$date)),]
env.data

frames <- lapply(dir("data/flow_cytometry/FCS/", full.names = TRUE), read.FCS)
names(frames) <- sapply(frames, keyword, "TUBE NAME")
fs <- as(frames, "flowSet")
sampleNames(fs)

phenoData(fs)

fsApply(fs, each_col, mean)

#f
# autoplot(transform(f, 
#                    "Pacific Blue-H" = asinh(`Pacific Blue-H`), 
#                    "APC-H" = asinh(`APC-H`)),
#          "Pacific Blue-H", "APC-H")

cell_counts <- data.frame(sample.id = names(frames), cell_count = NA)

# Controls
# 1) Just Hoechst
# 2) Just PY
# 3) Just Live-Dead efluor stain
# 6) Just Ecoli + beads
# 7) Ecoli + beads + PY + Hoechst

for(i in 1:length(fs)){
  f <- fs[[i]]
  f
  ft <- transform(f, 
                  "Pacific Blue-H" = asinh(`Pacific Blue-H`), 
                  "APC-H" = asinh(`APC-H`))
  presel <- list("Pacific Blue-H"=c(9,Inf), "APC-H"=c(8,10))
  
  # lg <- lymphGate(ft,
  #                 channels = c("Pacific Blue-H", "APC-H"), 
  #                 preselection = presel, scale = 100,
  #                 plot = F)
  
  rg <- rangeGate(ft, "Pacific Blue-H", alpha = "min",
                  refLine = 5, plot = TRUE)
  
  xyplot(`Pacific Blue-H`~`APC-H`, ft, filter = rg)
  
  fr <- flowCore::filter(ft,rg)
  count <- sum(fr@subSet)
  id <- fr@frameId
  cell_counts[i,] <- c(id, count)
}


cell_count_dat <- cell_counts[-c(1:7),] %>%
  mutate(cell_count = as.numeric(cell_count),
         sample.id = as.numeric(str_remove(sample.id, "UL_"))) %>% 
  mutate(cell_density = cell_count / 50000 * 10^6) # cells/50000 bead counts * 10^6 beads/ml


cell_plot <- env.data %>% left_join(cell_count_dat, by = "sample.id") %>% 
  filter(sample.id < 60) %>% 
  ggplot(aes(x = date, y = cell_density)) + 
  geom_point() +
  theme_minimal() +
  scale_y_log10(labels = label_scientific()) +
  scale_x_date(breaks = "3 month") + 
  labs(x = "",
       y = expression(paste("Cell density (cells ", ml^{-1},")")))

temp_plot <- env.data %>% left_join(cell_count_dat, by = "sample.id") %>% 
  filter(sample.id < 60) %>% 
  ggplot(aes(x = date, y = temp)) + 
  geom_point() +
  theme_minimal() +
  scale_x_date(breaks = "3 month") + 
  labs(x = "",
       y = expression(paste("Temperature (",~degree,"C)")))

cell_plot + temp_plot +
  plot_layout(nrow = 2, guides = "collect") +
  ggsave("figures/cell_density.pdf", width = 6, height = 3/4*6)
