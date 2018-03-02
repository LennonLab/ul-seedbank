require(codyn)
require(viridis)
require(betapart)


# never active
transients <- OTUs.REL.total[,which(colSums(OTUs.REL.active) == 0)]
# present, but not active
seedbank <- OTUs.REL.active == 0 & OTUs.REL.total != 0
seedbank <- seedbank * OTUs.REL.total
seedbank <- seedbank[,which(colSums(seedbank) > 0)]

seedbank.dist <- vegdist(seedbank, method = "euclidean")
seedbank.pcoa <- cmdscale(seedbank.dist, eig = T)
eigenvals(seedbank.pcoa)[1] / sum(eigenvals(seedbank.pcoa))
eigenvals(seedbank.pcoa)[2] / sum(eigenvals(seedbank.pcoa))
seedbank.traj <- seedbank.pcoa$points
colnames(seedbank.traj) <- c("Comp.1", "Comp.2")

png(filename = "figures/seedbank-dynamics.png", height = 8, width = 6, units = 'in', res = 300)
par(mfrow = c(3,1), mar = c(3,5,1,1), oma = c(2,2,2,2))
plot(env.data$sample.id, env.data$temp, type = 'l', ylab = "Temperature (C)", xaxt = "n", xlab = "", cex.lab = 1.3)
plot(rowSums(seedbank > 0), type = 'l', xaxt = "n", ylab = "Seedbank Richness", xlab = "", cex.lab = 1.3)
plot(seedbank.traj[,1], type = 'l', ylab = "Seedbank PC1", cex.lab = 1.3)
mtext("Time (Week)", side = 1, line = 3)
dev.off()

seedbank.traj.df <- cbind.data.frame(
  sample.id = design$sample.id[which(design$sample.type == "RNA")], seedbank.traj) 
pdf("figures/active_trajectory.pdf", bg = "white", height = 8, width = 10)
full_join(design[which(design$sample.type == "RNA"),], seedbank.traj.df) %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = sample.id)) + 
  geom_path(arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                          ends = "last", type = "closed")) + 
  scale_color_gradientn(colors = viridis(125), "Sample Number") + 
  geom_point(aes(x = seedbank.traj.df$Comp.1[1], y = seedbank.traj.df$Comp.2[1]),
             color = viridis(2)[1], size = 3) +
  geom_point(aes(x = seedbank.traj.df$Comp.1[123], y = seedbank.traj.df$Comp.2[123]),
             color = viridis(2)[2], size = 3) +
  theme_cowplot() + 
  xlab("PCoA 1") + 
  ylab("PCoA 2")
dev.off()


OTUs.dist <- vegdist(OTUs.REL, method = "euclidean")

# OTUs.pca <- princomp(OTUs.dist)
# OTUs.traj <- OTUs.pca$scores[,c(1,2)]

OTUs.traj <- cmdscale(OTUs.dist)
colnames(OTUs.traj) <- c("Comp.1", "Comp.2")
active.traj <- OTUs.traj[which(design$sample.type == "RNA"),]
total.traj <- OTUs.traj[which(design$sample.type == "DNA"),]

# active.traj.end <- active.traj[-1,]
# active.traj.start <- active.traj[-123,]
# active.traj.full <- cbind(active.traj.start, active.traj.end)
# colnames(active.traj.full) <- c("PC1.start", "PC2.start", "PC1.end", "PC2.end")

active.traj.df <- cbind.data.frame(
  sample.id = design$sample.id[which(design$sample.type == "RNA")], active.traj) 
pdf("figures/active_trajectory.pdf", bg = "white", height = 8, width = 10)
full_join(design[which(design$sample.type == "RNA"),], active.traj.df) %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = sample.id)) + 
  geom_path(arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                  ends = "last", type = "closed")) + 
  scale_color_gradientn(colors = viridis(125), "Sample Number") + 
  geom_point(aes(x = active.traj.df$Comp.1[1], y = active.traj.df$Comp.2[1]),
             color = viridis(2)[1], size = 3) +
  geom_point(aes(x = active.traj.df$Comp.1[123], y = active.traj.df$Comp.2[123]),
             color = viridis(2)[2], size = 3) +
  theme_cowplot() + 
  xlab("PCoA 1") + 
  ylab("PCoA 2")
dev.off()



total.traj.df <- cbind.data.frame(
  sample.id = design$sample.id[which(design$sample.type == "DNA")], total.traj) 
pdf("figures/total_trajectory.pdf", bg = "white", height = 8, width = 10)
full_join(design[which(design$sample.type == "DNA"),], total.traj.df) %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = sample.id)) + 
  geom_path(arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                          ends = "last", type = "closed")) + 
  scale_color_gradientn(colors = magma(125)) + 
  geom_point(aes(x = total.traj.df$Comp.1[1], y = total.traj.df$Comp.2[1]),
             color = magma(2)[1], size = 3) +
  geom_point(aes(x = total.traj.df$Comp.1[123], y = total.traj.df$Comp.2[123]),
             color = magma(2)[2], size = 3) +
  theme_cowplot()+
  xlab("PCoA 1")+
  ylab("PCoA 2")
dev.off()


matplot(cbind(total.traj[,1], active.traj[,1]), type = 'l', ylab = "PCoA1", xlab = "Week")
matplot(cbind(total.traj[,2], active.traj[,2]), type = 'l', ylab = "PCoA2", xlab = "Week")

active.dists <- as.matrix(OTUs.dist)[which(design$sample.type == "RNA"),which(design$sample.type == "RNA")]
total.dists <- as.matrix(OTUs.dist)[which(design$sample.type == "DNA"),which(design$sample.type == "DNA")]
plot(design$sample.id[which(design$sample.type == "RNA")], active.dists[,1])
plot(design$sample.id[which(design$sample.type == "DNA")], total.dists[,1])



## organize data structure for codyn
OTUs.long <- as.data.frame(OTUs.REL)
OTUs.long$mol <- design$sample.type
OTUs.long$number <- design$sample.id
OTUs.long <- gather(OTUs.long, key = otu, value = abundance, -mol, -number)

# Calculate turnover
turn.total <- turnover(df = OTUs.long, 
                 time.var = "number",
                 species.var = "otu",
                 abundance.var = "abundance",
                 replicate.var = "mol")
turn.appear <- turnover(df = OTUs.long, 
                       time.var = "number",
                       species.var = "otu",
                       abundance.var = "abundance",
                       replicate.var = "mol",
                       metric = "appearance")
turn.disappear <- turnover(df = OTUs.long, 
                       time.var = "number",
                       species.var = "otu",
                       abundance.var = "abundance",
                       replicate.var = "mol",
                       metric = "disappearance")
turn.total %>% mutate(Molecule = mol) %>%
  ggplot(aes(y = total, x = number, color = Molecule))+
  geom_line() + 
  ylab("Turnover") + 
  xlab("Time (Week)") +
  scale_color_brewer(type = "qual", palette = 6, direction = -1) +
  ggsave("figures/turnover-total.pdf", height = 6, width = 6, bg = "white")
ggplot(turn.appear, aes(y = appearance, x = number, color = mol))+
  geom_line()
ggplot(turn.disappear, aes(y = disappearance, x = number, color = mol))+
  geom_line()

dates <- left_join(design[-which(design$sample.id == 1),], env.data) %>% select(date)
turn.total$Season <- ifelse(month(dates$date) %in% c(5:10), "Summer", "Winter")
group_by(turn.total, mol, Season) %>% 
  summarise(
    mean = mean(total),
    cv = sd(total)/mean)

turn.total %>% group_by(mol) %>% summarize(avg.turn = mean(total))
# Mean Rank Shift
ul.mrs <- rank_shift(df = OTUs.long,
           time.var = "number",
           species.var = "otu",
           abundance.var = "abundance",
           replicate.var = "mol")
ul.mrs

ul.mrs$year <- as.numeric(colsplit(ul.mrs$year_pair, pattern = "-", names = c("year1", "year2"))[,2])
names(ul.mrs)[3] <- "Molecule"

# Create plot
pdf("figures/mean_rank_shift.pdf", bg = "white", height = 8, width = 8)
ggplot(ul.mrs, aes(x = year, y = MRS, color = Molecule)) + 
  geom_line(size = 1) +
  xlab("Time (Week)") + 
  ylab("Mean Rank Shift") +
  theme_cowplot() +
  scale_color_brewer(type = "qual", palette = 6, direction = -1)
dev.off()


# Does one plot type show higher or lower MRS, on average?
dates <- left_join(design[-which(design$sample.id == 1),], env.data) %>% select(date)

ul.mrs$Season <- ifelse(month(dates$date) %in% c(5:10), "Summer", "Winter")
group_by(ul.mrs, Molecule, Season) %>% 
  summarise(
    mean = mean(MRS),
    cv = sd(MRS)/mean)


# Rate change interval
ul.rci <- rate_change_interval(df = OTUs.long,
                     time.var = "number",
                     species.var = "otu", 
                     abundance.var = "abundance",
                     replicate.var = "mol")

ul.rci
rate.plot <- ggplot(ul.rci, aes(interval, distance)) +
  geom_point() + 
  facet_wrap(~mol) +
  theme(strip.text.x = element_text(size = 7)) +
  stat_smooth(method = "loess", se = F, size = 1) + 
  ylab("Hellinger Distance") + 
  xlab("Time Interval (Weeks)") +
  theme_minimal()
rate.plot


# Partition nestedness and turnover
OTUs.pa <- decostand(OTUs, method = "pa")
OTUs.betapart <- beta.pair(x = OTUs.pa, index.family = "sorensen")
OTUs.nest <- OTUs.betapart$beta.sne
OTUs.turn <- OTUs.betapart$beta.sim

OTUs.traj <- cmdscale(OTUs.turn)
colnames(OTUs.traj) <- c("Comp.1", "Comp.2")
active.traj <- OTUs.traj[which(design$sample.type == "RNA"),]
total.traj <- OTUs.traj[which(design$sample.type == "DNA"),]

# active.traj.end <- active.traj[-1,]
# active.traj.start <- active.traj[-123,]
# active.traj.full <- cbind(active.traj.start, active.traj.end)
# colnames(active.traj.full) <- c("PC1.start", "PC2.start", "PC1.end", "PC2.end")

active.traj.df <- cbind.data.frame(
  sample.id = design$sample.id[which(design$sample.type == "RNA")], active.traj) 
full_join(design[which(design$sample.type == "RNA"),], active.traj.df) %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = sample.id)) + 
  geom_path(arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                          ends = "last", type = "closed")) + 
  scale_color_gradientn(colors = viridis(2)) + 
  geom_point(aes(x = active.traj.df$Comp.1[1], y = active.traj.df$Comp.2[1]),
             color = viridis(2)[1], size = 3) +
  geom_point(aes(x = active.traj.df$Comp.1[123], y = active.traj.df$Comp.2[123]),
             color = viridis(2)[2], size = 3) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
  theme_minimal()



total.traj.df <- cbind.data.frame(
  sample.id = design$sample.id[which(design$sample.type == "DNA")], total.traj) 
full_join(design[which(design$sample.type == "DNA"),], total.traj.df) %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = sample.id)) + 
  geom_path(arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                          ends = "last", type = "closed")) + 
  scale_color_gradientn(colors = viridis(2)) + 
  geom_point(aes(x = total.traj.df$Comp.1[1], y = total.traj.df$Comp.2[1]),
             color = viridis(2)[1], size = 3) +
  geom_point(aes(x = total.traj.df$Comp.1[123], y = total.traj.df$Comp.2[123]),
             color = viridis(2)[2], size = 3) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
  theme_minimal()


matplot(cbind(total.traj[,1], active.traj[,1]), type = 'l', ylab = "PCoA1", xlab = "Week")
matplot(cbind(total.traj[,2], active.traj[,2]), type = 'l', ylab = "PCoA2", xlab = "Week")



## organize data structure for codyn
OTUs.long <- as.data.frame(seedbank)
OTUs.long$number <- design$sample.id[which(design$sample.type == "RNA")]
OTUs.long <- gather(OTUs.long, key = otu, value = abundance, -number)

# Calculate turnover
turn.total <- turnover(df = OTUs.long, 
                       time.var = "number",
                       species.var = "otu",
                       abundance.var = "abundance")
turn.appear <- turnover(df = OTUs.long, 
                        time.var = "number",
                        species.var = "otu",
                        abundance.var = "abundance",
                        metric = "appearance")
turn.disappear <- turnover(df = OTUs.long, 
                           time.var = "number",
                           species.var = "otu",
                           abundance.var = "abundance",
                           metric = "disappearance")

ggplot(turn.total, aes(y = total, x = number))+
  geom_line()
ggplot(turn.appear, aes(y = appearance, x = number))+
  geom_line()
ggplot(turn.disappear, aes(y = disappearance, x = number))+
  geom_line()

# Mean Rank Shift
ul.mrs <- rank_shift(df = OTUs.long,
                     time.var = "number",
                     species.var = "otu",
                     abundance.var = "abundance")
ul.mrs

ul.mrs$year <- as.numeric(colsplit(ul.mrs$year_pair, pattern = "-", names = c("year1", "year2"))[,2])
names(ul.mrs)[3] <- "Molecule"
# Create plot
pdf("figures/mean_rank_shift.pdf", bg = "white", height = 8, width = 8)
ggplot(ul.mrs, aes(x = year, y = MRS)) + 
  geom_line(size = 1) +
  xlab("Time (Week)") + 
  ylab("Mean Rank Shift") +
  theme_cowplot() +
  scale_color_brewer(type = "qual", palette = 6, direction = -1)
dev.off()
plot(rankshift.plot)

