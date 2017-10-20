par(mfrow = c(3,1))
matplot(OTUs.REL.active[,1], type = "l")
matplot(OTUs.REL.active[,2], type = "l")
matplot(OTUs.REL.active[,6], type = "l")

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
active.rich <- richness(OTUs.REL.active)
total.rich <- richness(OTUs.REL.sb)
matplot(cbind(active.rich, total.rich), type = 'l', col = c("red", "black"), lwd = 2, 
        ylab = "Richness", xlab = "Week")

active.alpha <- alpha(OTUs.REL.active)
total.alpha <- alpha(OTUs.REL.sb)
matplot(cbind(active.alpha, total.alpha), type = 'l', col = c("red", "black"), lwd = 2, 
        ylab = "D_alpha", xlab = "Week")


