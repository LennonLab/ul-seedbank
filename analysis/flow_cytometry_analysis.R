
library("flowCore")
library("flowViz")
library("flowStats")
library("ggcyto")
library("tidyverse")
library("patchwork")

design <- read.csv(file = "data/design.csv")

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

for(i in 1:length(fs)){
  f <- fs[[i]]
  ft <- transform(f, 
                  "Pacific Blue-H" = asinh(`Pacific Blue-H`), 
                  "APC-H" = asinh(`APC-H`))
  presel <- list("Pacific Blue-H"=c(9,Inf), "APC-H"=c(8,Inf))
  
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
         sample.id = as.numeric(str_remove(sample.id, "UL_")))


cell_plot <- env.dat %>% left_join(cell_count_dat, by = "sample.id") %>% 
  filter(sample.id < 60) %>% 
  ggplot(aes(x = sample.id, y = cell_count)) + 
  geom_point() +
  geom_line() +
  theme_minimal()

temp_plot <- env.dat %>% left_join(cell_count_dat, by = "sample.id") %>% 
  filter(sample.id < 60) %>% 
  ggplot(aes(x = sample.id, y = temp)) + 
  geom_point() +
  geom_line() +
  theme_minimal()

cell_plot + temp_plot +
  plot_layout(nrow = 2) 
