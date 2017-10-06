active.hel <- decostand(OTUs.active, method = "hellinger")
total.hel <- decostand(OTUs.total, method = "hellinger")

active.dist <- vegdist(active.hel, method = "euclidean")
total.dist <- vegdist(total.hel, method = "euclidean")

active.pca <- princomp(active.dist)
active.traj <- active.pca$scores[,c(1,2)]

total.pca <- princomp(total.dist)
total.traj <- total.pca$scores[,c(1,2)]

plot(active.traj, type = 'c', asp = 1, ylim = c(-1.5,1.5), xlim = c(-1.5, 1.5))
text(labels = seq(97,124,1), x = active.traj[,1], y = active.traj[,2])
points(total.traj, type = 'c', asp = 1, col = "red")
text(labels = seq(97,124,1), x = total.traj[,1], y = total.traj[,2], col = "red")
