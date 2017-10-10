require(codyn)

OTUs.hel <- decostand(OTUs, method = "hellinger")

OTUs.dist <- vegdist(OTUs.hel, method = "euclidean")

OTUs.pca <- princomp(OTUs.dist)
OTUs.traj <- OTUs.pca$scores[,c(1,2)]

active.traj <- OTUs.traj[which(design$sample.type == "RNA"),]
total.traj <- OTUs.traj[which(design$sample.type == "DNA"),]

plot(active.traj, type = 'c', asp = 1, ylim = c(-1.5,1.5), xlim = c(-1.5, 1.5))
text(labels = seq(97,124,1), x = active.traj[,1], y = active.traj[,2])
points(total.traj, type = 'c', asp = 1, col = "red")
text(labels = seq(97,124,1), x = total.traj[,1], y = total.traj[,2], col = "red")

matplot(cbind(total.traj[,1], active.traj[,1]))
active.dists <- as.matrix(OTUs.dist)[which(design$sample.type == "RNA"),which(design$sample.type == "RNA")]
total.dists <- as.matrix(OTUs.dist)[which(design$sample.type == "DNA"),which(design$sample.type == "DNA")]
plot(design$sample.number[which(design$sample.type == "RNA")], active.dists[,1])
plot(design$sample.number[which(design$sample.type == "RNA")], total.dists[,1])



## organize data structure for codyn
OTUs.long <- OTUs.REL
OTUs.long$mol <- design$sample.type
OTUs.long$number <- design$sample.number
OTUs.long <- gather(OTUs.long, key = otu, value = abundance, -mol, -number)

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

ggplot(turn.total, aes(y = total, x = number, color = mol))+
  geom_line()
ggplot(turn.appear, aes(y = appearance, x = number, color = mol))+
  geom_line()
ggplot(turn.disappear, aes(y = disappearance, x = number, color = mol))+
  geom_line()

# Mean Rank Shift
ul.mrs <- rank_shift(df = OTUs.long,
           time.var = "number",
           species.var = "otu",
           abundance.var = "abundance",
           replicate.var = "mol")
ul.mrs

ul.mrs$year <- abs(as.numeric(substr(ul.mrs$year_pair, 4, 7)))

# Create ggplot
rankshift.plot <- ggplot(ul.mrs, aes(x = year, y = MRS, color = mol)) + 
  geom_line(size = 1) + 
  xlab("Week") + 
  ylab("Mean Rank Shift") +
  theme_light()

plot(rankshift.plot)

# Does one plot type show higher or lower MRS, on average?
group_by(ul.mrs, mol) %>% 
  summarise(
    mean = mean(MRS),
    cv = sd(MRS)/mean)
