require(codyn)
require(viridis)

OTUs.dist <- vegdist(OTUs.REL, method = "euclidean")

OTUs.pca <- princomp(OTUs.dist)
OTUs.traj <- OTUs.pca$scores[,c(1,2)]

active.traj <- OTUs.traj[which(design$sample.type == "RNA"),]
total.traj <- OTUs.traj[which(design$sample.type == "DNA"),]

active.traj.end <- active.traj[-1,]
active.traj.start <- active.traj[-123,]
active.traj.full <- cbind(active.traj.start, active.traj.end)
colnames(active.traj.full) <- c("PC1.start", "PC2.start", "PC1.end", "PC2.end")

active.traj.df <- cbind.data.frame(
  sample.id = design$sample.id[which(design$sample.type == "RNA")], active.traj) 
full_join(design[which(design$sample.type == "RNA"),], active.traj.df) %>%
  ggplot(aes(x = Comp.1, y = Comp.2, color = sample.id)) + 
  geom_path(arrow = arrow(angle = 15, length = unit(0.1, "inches"),
                  ends = "last", type = "closed")) + 
  scale_color_gradientn(colors = viridis(2)) + 
  geom_point(aes(x = active.traj.df$Comp.1[1], y = active.traj.df$Comp.2[1]),
             color = viridis(2)[2], size = 3) +
  geom_point(aes(x = active.traj.df$Comp.1[122], y = active.traj.df$Comp.2[122]),
             color = viridis(2)[1], size = 3) +
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
  geom_point(aes(x = total.traj.df$Comp.1[122], y = total.traj.df$Comp.2[122]),
             color = viridis(2)[2], size = 3) +
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE) +
  theme_minimal()


  ggplot(aes(x = Comp.1, y = Comp.2, col = Sample.Num)) + 
  geom_line()

plot(active.traj, type = 'c', asp = 1, ylim = c(-2.5,2.5), xlim = c(-2.5, 2.5))
s = seq(1,122)
p.arrows(x1 = active.traj.full[s,1], y1 = active.traj.full[s,2],
         x2 = active.traj.full[s,3], y2 = active.traj.ƒull[s,4])
text(x = active.traj[,1], y = active.traj[,2])
points(total.traj, type = 'c', asp = 1, col = "red")
text(labels = seq(97,124,1), x = total.traj[,1], y = total.traj[,2], col = "red")

matplot(cbind(total.traj[,1], active.traj[,1]))
active.dists <- as.matrix(OTUs.dist)[which(design$sample.type == "RNA"),which(design$sample.type == "RNA")]
total.dists <- as.matrix(OTUs.dist)[which(design$sample.type == "DNA"),which(design$sample.type == "DNA")]
plot(design$sample.number[which(design$sample.type == "RNA")], active.dists[,1])
plot(design$sample.number[which(design$sample.type == "RNA")], total.dists[,1])



## organize data structure for codyn
OTUs.long <- as.data.frame(OTUs.REL)
OTUs.long$mol <- design$sample.type
OTUs.long$number <- design$sample.id
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

ul.mrs$year <- as.numeric(colsplit(ul.mrs$year_pair, pattern = "-", names = c("year1", "year2"))[,2])

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
