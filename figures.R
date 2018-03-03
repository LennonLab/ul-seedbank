pdf(file = "figures/temperature.pdf",width = 8, height = 8, bg = "white")
par(mar = c(5,5,4,3))
# Temperature
plot(data$sample.id, data$temp, xlim = c(0, 125), ylim = c(0, 33), type = "b", 
     pch = 22, bg = "grey", col = "black", cex = 1.2, ylab = "", xlab = "", 
     cex.lab = 1.5, las = 1, lwd = 2, yaxt = "n", xaxt = "n")
box(lwd=2)

mtext(expression(paste('Temp (',~degree~'C)', sep = '')), side = 2,
      outer = F, cex = 1.5, line = 3, adj = 0.5)  

# Major Axes
axis(side = 2, lwd.ticks = 2, tck = -0.02, cex.axis = 1.5, las = 1,
     labels = T, at = c(0, 10, 20, 30))
axis(side = 4, lwd.ticks = 2, cex.axis = 1, las = 1,
     labels = F, at = c(0, 10, 20, 30))
axis(side = 1, lwd.ticks = 2, tck = -0.02, cex.axis = 0.9, las = 1,
     labels = T, at = c(1, 25, 50, 75, 100, 125))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.5, las = 1,
     at = c(1, 25, 50, 75, 100, 125), labels = F)
mtext("Time (Week)", side = 1, cex = 2, 
      line = 3)
dev.off()



pdf(file = "figures/conductivity.pdf",width = 8, height = 8, bg = "white")
par(mar = c(5,5,4,3))

plot(data$sample.id, data$spc, xlim = c(0, 125), ylim = c(0.3, 0.7), type = "b", 
     pch = 22, bg = "grey", col = "black", cex = 1.2, ylab = "", xlab = "", 
     cex.lab = 1.5, las = 1, lwd = 2, yaxt = "n", xaxt = "n")
box(lwd=2)

mtext(expression(paste('Conductivity (',mu,'S/cm)', sep = '')), side = 2,
      outer = F, cex = 1.5, line = 3, adj = 0.5)  

# Major Axes
axis(side = 2, lwd.ticks = 2, tck = -0.02, cex.axis = 1.5, las = 1,
     labels = T, at = c(0.3, 0.5, 0.7))
axis(side = 4, lwd.ticks = 2, cex.axis = 1, las = 1,
     labels = F, at = c(0.3, 0.5, 0.7))
axis(side = 1, lwd.ticks = 2, tck = -0.02, cex.axis = 0.9, las = 1,
     labels = T, at = c(1, 25, 50, 75, 100, 125))
axis(side = 3, lwd.ticks = 2, cex.axis = 1.5, las = 1,
     at = c(1, 25, 50, 75, 100, 125), labels = F)
mtext("Time (Week)", side = 1, cex = 2, 
      line = 3)
dev.off()

