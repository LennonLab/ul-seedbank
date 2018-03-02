pdf("figures/taxa_dynamics.pdf", bg = "white", width = 8, height = 10)
par(mfrow = c(3,1), oma = c(2,2,2,2))
matplot(OTUs.REL.active[,1], type = "l", yaxt = "n", xaxt = "n",
        ylab = "", xlab = "", lwd = 2)
box(lwd=2)
mtext(expression(italic("Anabaena sp.")), side = 2,
      outer = F, cex = 1.5, line = 3, adj = 0.5)  

# Major Axes
axis(side = 2, lwd.ticks = 2, tck = -0.02, cex.axis = 1.5, las = 1,
     labels = T)
axis(side = 4, lwd.ticks = 2, cex.axis = 1, las = 1,
     labels = F)
axis(side = 1, lwd.ticks = 2, tck = -0.02, cex.axis = 1, las = 1,
     labels = F, at = c(1, 25, 50, 75, 100, 125))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.5, las = 1,
     at = c(1, 25, 50, 75, 100, 125), labels = F)

matplot(OTUs.REL.active[,2], type = "l", yaxt = "n", xaxt = "n",
        ylab = "", xlab = "", lwd = 2)
box(lwd=2)
mtext(expression(italic("Limnohabitans sp.")), side = 2,
      outer = F, cex = 1.5, line = 3, adj = 0.5)  
# Major Axes
axis(side = 2, lwd.ticks = 2, tck = -0.02, cex.axis = 1.5, las = 1,
     labels = T)
axis(side = 4, lwd.ticks = 2, cex.axis = 1, las = 1,
     labels = F)
axis(side = 1, lwd.ticks = 2, tck = -0.02, cex.axis = 1, las = 1,
     labels = F, at = c(1, 25, 50, 75, 100, 125))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.5, las = 1,
     at = c(1, 25, 50, 75, 100, 125), labels = F)


matplot(OTUs.REL.active[,6], type = "l", yaxt = "n", xaxt = "n",
        ylab = "", xlab = "", lwd = 2)
box(lwd=2)
mtext(expression(italic("Terrimicrobium sp.")), side = 2,
      outer = F, cex = 1.5, line = 3, adj = 0.5)  

# Major Axes
axis(side = 2, lwd.ticks = 2, tck = -0.02, cex.axis = 1.5, las = 1,
     labels = T)
axis(side = 4, lwd.ticks = 2, cex.axis = 1, las = 1,
     labels = F)
axis(side = 1, lwd.ticks = 2, tck = -0.02, cex.axis = 1.5, las = 1,
     labels = T, at = c(1, 25, 50, 75, 100, 125))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.5, las = 1,
     at = c(1, 25, 50, 75, 100, 125), labels = F)
mtext("Time (Weeks)", side = 1, line = 3.5, cex = 1.5)
dev.off()

par(mfrow = c(3,1))
matplot(OTUs.REL.total[,1], type = "l")
matplot(OTUs.REL.total[,2], type = "l")
matplot(OTUs.REL.total[,6], type = "l")

require(zoo)
temp.ts <- zoo(x = env.data$temp, order.by = env.data$date)
otu.ts <- zoo(x = OTUs.REL.active, order.by = env.data$date)
head(otu.ts)
plot(temp.ts)
otu.ts[,1]

ccf(temp.ts, otu.ts[,1], na.action = na.pass)
require(WaveletComp)

cbind.data.frame(env.data)
temp.wave <- analyze.wavelet(na.omit(env.data[,1:3]), my.series = 3, loess.span = 0)
temp.wave
wt.image(temp.wave)
wt.avg(temp.wave, siglvl = 0.05, sigcol = "red")
plot(temp.wave$Period, temp.wave$Power.avg, type = 'l')

require(hydrostats)
Colwells(temp.ts)
Colwells(ts.format(as.data.frame(env.data[,c(2:3)]), format = "%Y-%m-%d", cols = c(1, 2)))

dates <- env.data[which(env.data$sample.id %in% design$sample.id),2]
sp.seasonality <- function(taxon.abund, dates){
  M <- Colwells(
    ts.format(
      data.frame(dates = dates, Q = taxon.abund), 
      format = "%Y-%m-%d",
      cols = c(1, 2)),
    boundaries = "log_class_size",
    indices.only = T)
  return(M$M)
}
sp.seasonality(taxon.abund = OTUs.REL.active[,1], dates = dates)


richness <- function(M){
  rowSums(decostand(M, method = "pa"))
}

alpha <- function(M){
  apply(M, MARGIN = 1, d)
}

div <- function(M){
  diversity(M, index = "invsimpson")
}

active.rich <- richness(OTUs.REL.active)
total.rich <- richness(OTUs.REL.total)
matplot(cbind(active.rich, total.rich), type = 'l', col = c("red", "black"), lwd = 2, 
        ylab = "Richness", xlab = "Week")
cbind.data.frame(RNA = active.rich, DNA = total.rich) %>% 
  rowid_to_column(var = "sample.id") %>% left_join(env.data) %>% select(RNA, DNA, date) %>%
  gather(Molecule, Richness, -date) %>%
  ggplot(aes(x = as.Date(date), y = Richness, col = Molecule)) + 
  geom_line() +
  xlab("Date") +
  theme(axis.text.x=element_text(angle=60, hjust=1)) + 
  scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y") + 
  scale_color_brewer(type = "qual", palette = 6, direction = -1) +
  ggsave("figures/richness.pdf", height = 6, width = 6, bg = "white")

active.alpha <- alpha(OTUs.REL.active)
total.alpha <- alpha(OTUs.REL.sb)
matplot(cbind(active.alpha, total.alpha), type = 'l', col = c("red", "black"), lwd = 2, 
        ylab = "D_alpha", xlab = "Week")


weatherdat <- read.csv(file = "data/noaa-climate.csv")
weatherdat %>% filter(NAME == "BLOOMINGTON INDIANA UNIVERSITY, IN US") %>%
  mutate(Date = as_date(DATE)) %>% mutate(Week = weeks(Date)) %>%
  select(Date, PRCP) %>% 
  ggplot(aes(y = PRCP, x = Date)) + 
  geom_line()


weatherdat %>% filter(NAME == "BLOOMINGTON INDIANA UNIVERSITY, IN US") %>%
  mutate(Date = as_date(DATE)) %>% select(Date, PRCP) %>% 
  ggplot(aes(y = PRCP, x = Date)) + 
  geom_line()
samplingdates <- as_date(env.data$date)
weatherdates <- weatherdates$Date
weatherdat.subset <- weatherdat[which(weatherdates %in% samplingdates),]
weatherdat.subset %>% select(Date,)
