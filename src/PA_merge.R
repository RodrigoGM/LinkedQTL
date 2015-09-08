## Analysis of cf2_o
## Working directories
setwd("./analysis/")

## libraries
library(RColorBrewer)
library(gplots)
library(lattice)

## merging real with estimated data
##  To merge only based on analysis and on qtl model (q1, q2)
## run ../src/MergeOutputs.sh or run from within R
system("bash ../src/MergeOutputs.sh")

mqm.q1 <- read.csv("MQM_q1_summary.csv")
mqm.q2 <- read.csv("MQM_q2_summary.csv")

im.q1 <-  read.csv("S1_q1_summary.csv")
im.q2 <-  read.csv("S1_q2_summary.csv")

modelq1 <- read.csv("../cf2_d/CF2-q1_Model.csv")
modelq2 <- read.csv("../cf2_d/CF2-q2_Model.csv")

colnames(modelq1)[1] <- "cross.no"
colnames(modelq2)[1] <- "cross.no"
modelq1$model <- 1
modelq2$model <- 2

mqm.q1 <- merge(mqm.q1, modelq1, by = "cross.no", all.x = TRUE)
mqm.q2 <- merge(mqm.q2, modelq2, by = "cross.no", all.x = TRUE)

im.q1 <- merge(im.q1, modelq1, by = "cross.no", all.x = TRUE)
im.q2 <- merge(im.q2, modelq2, by = "cross.no", all.x = TRUE)

mqm <- rbind(mqm.q1, mqm.q2)
im <- rbind(im.q1, im.q2)


## adding descriptive columns
mqm$corr[mqm$phenotype == "phenotype.1"] <- 1
mqm$corr[mqm$phenotype == "phenotype.2"] <- .90
mqm$corr[mqm$phenotype == "phenotype.3"] <- .60
mqm$corr[mqm$phenotype == "phenotype.4"] <- .30
mqm$corr[mqm$phenotype == "phenotype.5"] <- 0

im$corr[im$phenotype == "phenotype.1"] <- 1
im$corr[im$phenotype == "phenotype.2"] <- .90
im$corr[im$phenotype == "phenotype.3"] <- .60
im$corr[im$phenotype == "phenotype.4"] <- .30
im$corr[im$phenotype == "phenotype.5"] <- 0

## estimating distance from real qtl
mqm$dist <- sapply(1:length(mqm[,1]), function(i) min(abs(as.numeric(mqm$peak[i]) - as.numeric(mqm$loc1[i])), abs(as.numeric(mqm$peak[i]) - as.numeric(mqm$loc2[i]))))
im$dist <- sapply(1:length(im[,1]), function(i) min(abs(as.numeric(im$peak[i]) - as.numeric(im$loc1[i])), abs(as.numeric(im$peak[i]) - as.numeric(im$loc2[i]))))

## estimating average location of qtl.  For plotting based on a single coordinate
mqm$loc <- sapply(1:length(mqm[,1]), function(i) mean(c(mqm$loc1[i], mqm$loc2[i])))
im$loc <-  sapply(1:length(im[,1]), function(i) mean(c(im$loc1[i], im$loc2[i])))

## estimating confidence interval length
mqm$ci.length <- as.numeric(mqm$uci) - as.numeric(mqm$lci)
im$ci.length <- as.numeric(im$uci) - as.numeric(im$lci)

## creating nice pretty labels
mqm$model.t <- ifelse(mqm$model == 1, "a1 = 7.5 ; a2 = 7.5", "a1 = 5 ; a2 = 10")
mqm$corr.t <- paste("r =", mqm$corr)

im$model.t <- ifelse(im$model == 1, "a1 = 7.5 ; a2 = 7.5", "a1 = 5 ; a2 = 10")
im$corr.t <- paste("r =", im$corr)

mqm$nind.t <- paste("N =", mqm$nind)
im$nind.t <- paste("N =", im$nind)

## first generation plots
pdf("../figures/first_plots.pdf")
plotmeans(dist ~ phenotype, data = mqm)
plotmeans(ci.length ~ phenotype, data = mqm)
dev.off()

## writing out main merged output files
write.table(mqm, file = "MQM_summary.tab", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(im, file = "IM_summary.tab", sep = "\t", quote = FALSE, row.names = FALSE)

qtls <- data.frame(cross.no = matrix(NA, nrow = nrow(mqm) + nrow(im)))

qtls$cross.no <- c(mqm$cross.no, im$cross.no)
qtls$name <- c(mqm$name, im$name)
qtls$phenotype <- c(mqm$phenotype, im$phenotype)

qtls <- rbind(mqm[, intersect(names(mqm), names(im))], im[, intersect(names(mqm), names(im))])
qtls$lod <- c(mqm$lod, im$peak.lod)
qtls$known.dist <- c(mqm$k.dist, im$known.dist)
qtls$pLOD <- c(mqm$pLOD, rep(NA, nrow(im)))
qtls$n.gen <- c(mqm$n.gen, rep(NA, nrow(im)))
qtls$n.qtl <- c(mqm$n.qtl, rep(NA, nrow(im)))
qtls$peak.lod <- c(mqm$lod, im$peak.lod)
qtls$lci.lod <- c(rep(NA, nrow(mqm)), im$lci.lod)
qtls$uci.lod <- c(rep(NA, nrow(mqm)), im$uci.lod)

write.table(qtls[1:564578,], file = "QTL_summary1.tab", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(qtls[564579:1129155,], file = "QTL_summary2.tab", sep = "\t", quote = FALSE, row.names = FALSE)

save.image("LinkedQTL.rda")


qtls$method <- as.character(qtls$method)
qtls$method<- gsub("stq", "QMS", qtls$method)
qtls$detected <- sapply(1:length(qtls[,1]), function(i) ifelse(abs(qtls$peak[i] - qtls$loc1[i]) < abs(qtls$peak[i] - qtls$loc2[i]), "a1", "a2"))

## Reordering...
qtls$nind.t <- factor(qtls$nind.t, levels = c("N = 1500", "N = 1000", "N = 500", "N = 250"))
qtls$model.t <- factor(qtls$model.t, levels = c("a1 = 7.5 ; a2 = 7.5", "a1 = 5 ; a2 = 10"))
qtls$method <- factor(qtls$method, levels = c("IM", "CIM", "QMS"))
qtls$model.ta <- factor(qtls$model.t, levels = c("a1 = 5 ; a2 = 10", "a1 = 7.5 ; a2 = 7.5"))

qtls$cd <- NA
###qtls$cd[qtls$detected == "a1" && qtls$lci <= qtls$loc1 && qtls$loc1 <= qtls$uci] <-  1
###qtls$cd[qtls$detected == "a2" && qtls$lci <= qtls$loc2 && qtls$loc2 <= qtls$uci] <- 1
qtls$fd <- NA
###qtls$fd[qtls$cd == 1] <- 0
###qtls$fd[qtls$cd == 0] <- 1
## check for cd and fd in excel nd re read file.
#write.table(qtls, file = "QTL_summary.tab", sep = "\t", quote = FALSE, row.names = FALSE)




