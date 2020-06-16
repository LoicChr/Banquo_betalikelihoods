obs.comm <- read.csv("data/kettlehole_dataset/community.cover.15sp.csv", header=T, row.names = 1)
load("data/kettlehole_dataset/Data_mod.Rdata")

data <- KH.data
data$QuadID <- droplevels(data$QuadID)
data$poros <- log(data$poros)

t.avg <- ddply(KH.data, c("species"), summarise,
               ht = mean(height),maxht  = quantile(height,0.975), sla = mean(sla, na.rm = T))
t.avg$maxht <- log(t.avg$maxht)
t.avg$sla <- log(t.avg$sla)

dudi.tr <- dudi.pca(t.avg[,c("maxht", "sla")], nf = 2, scannf = F)
t.avg$Axis1 <- dudi.tr$li[,1]
t.avg$Axis2 <- dudi.tr$li[,2]

# Load traitspace object
load("data/tra2003/traitspace_unnorm.Rdata")
P_S_E_tr <- as.matrix(P_S_E_tra)

comm.red <- obs.comm[row.names(P_S_E_tr), colnames(P_S_E_tr)]
# comm.red <- as.matrix(sweep(comm.red, 1, rowSums(comm.red), "/")) +0.01
# comm.red <- as.matrix(sweep(comm.red, 1, rowSums(comm.red), "/")) 

load("data/ozab_data.Rdata")
