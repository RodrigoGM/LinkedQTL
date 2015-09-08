## Analysis of cf2_o
## Working directories
setwd("./analysis/")

source("../src/PA_merge.R")

cols <- brewer.pal(9, "Set1")
colb <- brewer.pal(3, "Dark2")

out.d <- glm(dist ~ corr + model + nind + nind:model + method + method:model +
  method:nind + method:nind:model + loc + loc:model + loc:nind:model +
  loc:method:model, data = qtls)

pdf("glm diagnostics dist.pdf", width = 11, height = 7)
plot(out.d, pch = 19, cex = .5)
dev.off()


out.l <- glm(ci.length ~ model + nind + nind:model + nmar + nind:nmar +
  nind:nmar:model + method + corr:method + method:model + method:nind +
  method:nind:model + method:nind:nmar + loc + loc:model + loc:nind +
  loc:nind:model + loc:nind:nmar + loc:method + corr:loc:method +
  loc:method:model + loc:method:nind, data = qtls)

pdf("glm diagnostics ci_length.pdf", width = 11, height = 7)
plot(out.l, pch = 19, cex = .5)
dev.off()


pdf("dist and ci length, Bopxplots.pdf", width = 8, height = 12)
par(mfrow = c(3, 1), bty = "l")
boxplot(dist ~ model.t*nind.t, data = qtls[qtls$method == "QMS",], col = cols[2],
  main = "Distance with MQM", cex = 2, pch = ".")
boxplot(dist ~ model.t*nind.t, data = qtls[qtls$method == "QMS",],
  col = cols[3], main = "Distance with IM", cex = 2, pch = ".")
boxplot(dist ~ model.t*nind.t, data = qtls[qtls$method == "QMS",],
  col = cols[4], main = "Distance with CIM", cex = 2, pch = ".")

boxplot(ci.length ~ model.t*nind.t, data = qtls[qtls$method == "QMS",], col = cols[2],
  main = "CI Length with MQM", cex = 2, pch = ".")
boxplot(ci.length ~ model.t*nind.t, data = qtls[qtls$method == "QMS",],
  col = cols[3], main = "CI Length with IM", cex = 2, pch = ".")
boxplot(ci.length ~ model.t*nind.t, data = qtls[qtls$method == "QMS",],
  col = cols[4], main = "CI Length with CIM", cex = 2, pch = ".")

dev.off()

##Histogram of QTL locations
tiff("Histogram of QTL Locations.tiff", width = 1200, height = 500, units = "px",
  pointsize = 14)
par(mar = c(4.2, 4, 2, 1)+.1, mfrow = c(2, 2), cex.lab = 1.3, cex.axis = 1.2)
hist(modelq1$loc1, breaks = 500, xlim = c(0,100), main = "", xlab = "")
hist(modelq2$loc1, breaks = 500, xlim = c(0,100), main = "", xlab = "", ylab = "")
hist(modelq1$loc2, breaks = 500, xlim = c(0,100), main = "", xlab = "QTL Location (cM)")
hist(modelq2$loc2, breaks = 500, xlim = c(0,100), main = "", xlab = "QTL Location (cM)", ylab = "")
dev.off()


             
#pdf("Dist and CI length PlotMeans.pdf", width = 10, height =10)

tiff("Dist vs nind or nmar by method.tiff", width = 1600, height = 1000, units = "px",
  pointsize = 15)
  
nf <- layout(mat = matrix(c(4, 3, 1, 2), nrow = 2, ncol = 2), widths = 1, heights = 1)
layout.show(nf)
par(mar = c(4.4, 5, 2, 1.3)+.1, lwd = 3, cex.axis = 2, cex.lab = 2.1, bty = "l")
ylimd = c(0, 5)
ylimci = c(0, 10)
# 1 #
plotmedians(dist ~ nind, data = qtls[qtls$model == 1 & qtls$method == "QMS",], 
  col = cols[3], pch = 19, n.label = FALSE, ylim = ylimd, cex = 2.3,
  xlab = "Cross Size (N)", ylab = "Distance (N)",
  cex.lab = 2, cex.axis = 2)
plotmedians(dist ~ nind, data = qtls[qtls$model == 1 & qtls$method == "IM",], 
  col = cols[1], pch = 19, cex = 2.3, n.label = FALSE, add = TRUE, xaxt = "n", yaxt = "n")
plotmedians(dist ~ nind, data = qtls[qtls$model == 1 & qtls$method == "CIM",], 
  col = cols[2], pch = 19, cex = 2.3, n.label = FALSE, add = TRUE, xaxt = "n", yaxt = "n")

plotmedians(dist ~ nind, data = qtls[qtls$model == 2 & qtls$method == "QMS",], 
  col = cols[3], pch = 18, n.label = FALSE, lty = 2, cex = 2.3, 
  ylim = c(0, 12.5), add = TRUE,
  xlab = "", ylab = "")
plotmedians(dist ~ nind, data = qtls[qtls$model == 2 & qtls$method == "IM",], 
  col = cols[1], pch = 18, n.label = FALSE, lty = 2, cex = 2.3, add = TRUE, 
  xaxt = "n", yaxt = "n")
plotmedians(dist ~ nind, data = qtls[qtls$model == 2 & qtls$method == "CIM",], 
  col = cols[2], pch = 18, n.label = FALSE, lty = 2, cex = 2.3, add = TRUE, 
  xaxt = "n", yaxt = "n")

#legend("topright", legend = c("a1 = 5 | a2 = 5 | QMS", "a1 = 5 | a2 = 5 | IM", "a1 = 5 | a2 = 5 | CIM",
#                              "a1 = 5 | a2 = 10 | QMS", "a1 = 5 | a2 = 10 | IM", "a1 = 5 | a2 = 10 | CIM"),
#        lty = c(1,1,1, 2,2,2), pch = c(19, 19, 19, 18, 18, 18), col = rep(cols[1:3],2),
#        cex = 1.3)

# 2 #
plotmedians(ci.length ~ nind, data = qtls[qtls$model == 1 & qtls$method == "QMS",], 
  col = cols[3], pch = 19, n.label = FALSE, ylim = ylimci, cex = 2.3,  
  xlab = "Cross Size (N)", ylab = "CI Length (cM)")
plotmedians(ci.length ~ nind, data = qtls[qtls$model == 1 & qtls$method == "IM",], 
  col = cols[1], pch = 19, cex = 2.3, 
  n.label = FALSE, add = TRUE, xaxt = "n", yaxt = "n")
plotmedians(ci.length ~ nind, data = qtls[qtls$model == 1 & qtls$method == "CIM",], 
  col = cols[2], pch = 19, cex = 2.3, 
  n.label = FALSE, add = TRUE, xaxt = "n", yaxt = "n")

plotmedians(ci.length ~ nind, data = qtls[qtls$model == 2 & qtls$method == "QMS",], 
  col = cols[3], pch = 18, cex = 2.3, 
  n.label = FALSE, lty = 2, ylim = c(0, 12.5), add = TRUE, xaxt = "n", yaxt = "n",
  xlab = "", ylab = "")
plotmedians(ci.length ~ nind, data = qtls[qtls$model == 2 & qtls$method == "IM",], 
  col = cols[1], pch = 18, cex = 2.3, 
  n.label = FALSE, lty = 2, add = TRUE, xaxt = "n", yaxt = "n")
plotmedians(ci.length ~ nind, data = qtls[qtls$model == 2 & qtls$method == "CIM",], 
  col = cols[2], pch = 18, cex = 2.3, 
  n.label = FALSE, lty = 2, add = TRUE, xaxt = "n", yaxt = "n")

legend("topright", legend = c("a1 = 5 | a2 = 5 | IM", "a1 = 5 | a2 = 5 | CIM", "a1 = 5 | a2 = 5 | QMS", 
                              "a1 = 5 | a2 = 10 | IM", "a1 = 5 | a2 = 10 | CIM", "a1 = 5 | a2 = 10 | QMS"),
        lty = c(1,1,1, 2,2,2), pch = c(19, 19, 19, 18, 18, 18), col = rep(cols[c(1,2,3)],2),
        cex = 2)

#dev.off()

# 3 #
#pdf("CI Length vs Markers plotmedians.pdf", width = 8, height = 8)
#par(lwd = 2, cex.axis = 1.3, cex.lab = 1.3, bty = "l")
plotmedians(ci.length ~ nmar, data = qtls[qtls$method == "QMS" & qtls$model == 1,], 
  col = cols[3], ylim = ylimci, #main = "CI Length vs Marker Density",
  pch = 19, cex = 2.3, 
  xlab = "Marker Density", ylab = "CI Length (cM)", n.label = FALSE)
plotmedians(ci.length ~ nmar, data = qtls[qtls$method == "IM" & qtls$model == 1,],  
  pch = 19, col = cols[1], cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")
plotmedians(ci.length ~ nmar, data = qtls[qtls$method == "CIM" & qtls$model == 1,],
  pch = 19, col = cols[2], cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")

plotmedians(ci.length ~ nmar, data = qtls[qtls$method == "QMS" & qtls$model == 2,], 
  pch = 18, col = cols[3], lty = 2, cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n",
  ylab = "", xlab = "")
plotmedians(ci.length ~ nmar, data = qtls[qtls$method == "IM" & qtls$model == 2,],  
  pch = 18, col = cols[1], lty = 2, cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")
plotmedians(ci.length ~ nmar, data = qtls[qtls$method == "CIM" & qtls$model == 2,],
  pch = 18, col = cols[2], lty = 2, cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")

#legend("topright", legend = c("a1 = 5 | a2 = 5 | QMS", "a1 = 5 | a2 = 5 | IM", "a1 = 5 | a2 = 5 | CIM",
#                              "a1 = 5 | a2 = 10 | QMS", "a1 = 5 | a2 = 10 | IM", "a1 = 5 | a2 = 10 | CIM"),
#  col = cols[1:3], pch = c(rep(19, 3), rep(18, 3)), lty = c(rep(1, 3), rep(2,3)), 
#    lwd = 2, bty = "o", cex = 1.2)
  

# 4 #
plotmedians(dist ~ nmar, data = qtls[qtls$method == "QMS" & qtls$model == 1,], 
  col = cols[3], ylim = ylimd, #main = "CI Length vs Marker Density",
  pch = 19, cex = 2.3, 
  xlab = "Number of Markers (Md)", ylab = "Distance (cM)", n.label = FALSE)
plotmedians(dist ~ nmar, data = qtls[qtls$method == "IM" & qtls$model == 1,],  
  pch = 19, col = cols[1], cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")
plotmedians(dist ~ nmar, data = qtls[qtls$method == "CIM" & qtls$model == 1,],
  pch = 19, col = cols[2], cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")

plotmedians(dist ~ nmar, data = qtls[qtls$method == "QMS" & qtls$model == 2,], 
  pch = 18, col = cols[3], lty = 2, cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n") #, ylim = c(4, 10) )
plotmedians(dist ~ nmar, data = qtls[qtls$method == "IM" & qtls$model == 2,],  
  pch = 18, col = cols[1], lty = 2, cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")
plotmedians(dist ~ nmar, data = qtls[qtls$method == "CIM" & qtls$model == 2,],
  pch = 18, col = cols[2], lty = 2, cex = 2.3, 
  add = TRUE, n.label = FALSE, xaxt = "n", yaxt = "n")

#legend("topright", legend = c("a1 = 5 | a2 = 5 | QMS", "a1 = 5 | a2 = 5 | IM", "a1 = 5 | a2 = 5 | CIM",
#                              "a1 = 5 | a2 = 10 | QMS", "a1 = 5 | a2 = 10 | IM", "a1 = 5 | a2 = 10 | CIM"),
#  col = cols[1:3], pch = c(rep(19, 3), rep(18, 3)), lty = c(rep(1, 3), rep(2,3)), 
#    lwd = 2, bty = "o", cex = 1.2)


dev.off()


### LATICE PLOTS TO UNDERSTAND LOCATION EFFECTS

nind.key = list(text = list(c("N = 250", "N = 500", "N = 1000", "N = 1500"), font = 2, cex = 1.3),
            points = list(col = cols[1:4],pch = 19, cex = 1.5), space = "right"
          )

mar.key = list(text = list(c("Md = 96", "Md = 384", "Md = 1152", "Md = 1536"), font = 2, cex = 1.3),
            points = list(col = cols[1:4], pch = 19, cex = 1.5), space = "right"
          )

wi = 860
he = 500
xlim1 = c(5, 95)
ylim1 = c(0, 100)
##pdf("Peak, dist, ci length vs loc.pdf", width = 11, height = 8)
# Peak
tiff("Peak vs loc All.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)

xyplot(peak ~ loc|method*model.ta, data = qtls,
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nind, col = cols, main = "", pch = ".",
  xlim = xlim1, ylim = ylim1,
  cex = 2, cex.axis = 1.6, cex.lab = 2,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
#  auto.key = TRUE)
  key = nind.key)
dev.off()

tiff("Peak vs loc by nind QMS.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)

xyplot(peak ~ loc|model.t*nind.t, data = qtls[qtls$method == "QMS",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "",
  pch = ".", cex = 2, cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()

tiff("Peak vs loc by nind IM.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(peak ~ loc|model.t*nind.t, data = qtls[qtls$method == "IM",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "",
  pch = ".", cex = 2, cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()

tiff("Peak vs loc by nind CIM.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(peak ~ loc|model.t*nind.t, data = qtls[qtls$method == "CIM",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "",
  pch = ".", cex = 2, cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
#  auto.key = TRUE)
  key = mar.key)

dev.off()


# Confidence Interval Length

tiff("CI Length vs loc All.tiffa", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(ci.length ~ loc|method*model.ta, data = qtls,
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nind, col = cols, main = "", pch = ".",cex = 2, 
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Length of Confidence Interval (cM)",
#  auto.key = TRUE)
  key = nind.key)
dev.off()

tiff("CI Length vs loc by nind QMS.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(ci.length ~ loc|model.t*nind.t, data = qtls[qtls$method == "QMS",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "", pch = ".",cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Length of Confidence Interval (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()

tiff("CI Length vs loc by nind IM.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(ci.length ~ loc|model.t*nind.t, data = qtls[qtls$method == "IM",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "", pch = ".", cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Length of Confidence Interval (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()

tiff("CI Length vs loc by nind CIM.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(ci.length ~ loc|model.t*nind.t, data = qtls[qtls$method == "CIM",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "", pch = ".", cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Length of Confidence Interval (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()


# Accuracy

tiff("Dist vs loc All.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(dist ~ loc|method*model.ta, data = qtls,
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nind, col = cols, main = "", pch = ".",cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Distance from Real QTL Location (cM)",
#  auto.key = TRUE)
  key = nind.key)
dev.off()

tiff("Dist vs loc by nind QMS.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(dist ~ loc|model.t*nind.t, data = qtls[qtls$method == "QMS",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "", pch = ".", cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Distance from Real QTL Location (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()

tiff("Dist vs loc by nind IM.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(dist ~ loc|model.t*nind.t, data = qtls[qtls$method == "IM",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "", pch = ".", cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Distance from Real QTL Location (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()

tiff("Dist vs loc by nind CIM.tiff", width = wi, height = he, units = "px", 
    pointsize = 20)
xyplot(dist ~ loc|model.t*nind.t, data = qtls[qtls$method == "CIM",],
  strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE), 
  groups = nmar, col = cols, main = "", pch = ".", cex = 2,
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Distance from Real QTL Location (cM)",
#  auto.key = TRUE)
  key = mar.key)
dev.off()


## Peaks by QTL identified correctly

tiff("CD Peak vs loc by detected QMS.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(peak ~ loc|model.t*nind.t, qtls[qtls$method == "QMS" & qtls$cd == 1,], 
  groups = detected, col = cols[3:4],
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("a1", "a2"), font = 2, cex = 1.5),
            points = list(col = cols[3:4], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()
  
tiff("CD Peak vs loc by detected IM.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(peak ~ loc|model.t*nind.t, qtls[qtls$method == "IM" & qtls$cd == 1,], 
  groups = detected, col = cols[3:4],
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("a1", "a2"), font = 2, cex = 1.5),
            points = list(col = cols[3:4], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()

tiff("CD Peak vs loc by detected CIM.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
#layout(mat = matrix(c(1,2), nrow = 1, ncol = 2))
xyplot(peak ~ loc|model.t*nind.t, qtls[qtls$method == "CIM" & qtls$cd == 1,], 
  groups = detected, col = cols[3:4],
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("a1", "a2"), font = 2, cex = 1.5),
            points = list(col = cols[3:4], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()



## By correct detection or incorrect detection
xyplot(peak ~ loc|method*model.ta, qtls, 
  groups = cd, col = cols,
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )



tiff("Peak vs loc by CD QMS.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(peak ~ loc|model.t*nind.t, qtls[qtls$method == "QMS",], 
  groups = cd, col = cols,
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()
  
tiff("Peak vs loc by CD IM.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(peak ~ loc|model.t*nind.t, qtls[qtls$method == "IM",], 
  groups = cd, col = cols,
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()

tiff("Peak vs loc by CD CIM.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
#layout(mat = matrix(c(1,2), nrow = 1, ncol = 2))
xyplot(peak ~ loc|model.t*nind.t, qtls[qtls$method == "CIM",], 
 strip = strip.custom(par.strip.text = list(cex = 1.0, font = 2, col = "black"),
    strip.names=FALSE, strip.levels=TRUE),  groups = cd, col = cols, 
  pch = ".", cex = 2, #auto.key = TRUE
  cex.axis = 1.6, cex.lab = 2,
  xlim = xlim1, ylim = ylim1,
  xlab = "Average Real Location of QTL (cM)", ylab = "Estimated Location of QTL (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()



tiff("Distribution of Distance by nind.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
histogram(~dist|method*model.ta, data = qtls, 
  groups = nind.t, breaks = 100, 
  col = cols[4:1], xlim = c(0,25), ylim = c(0, 100),
  xlab = "cM",
  panel = function(...) {
  panel.grid(h = 4, v = 4)
  panel.superpose(...)
  },    panel.groups = function(...) {
    panel.histogram(..., border = col)  
    },
  key = list(text = nind.key$text, space = "right", 
    rectangles = list(col = cols[1:4]))
  )

dev.off()

tiff("Distribution of CI length by nind.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
histogram(~ci.length|method*model.ta, data = qtls, 
  groups = nind.t, breaks = 100, 
  col = cols[4:1], xlim = c(0,25), ylim = c(0, 50), 
  xlab = "cM",
  panel = function(...) {
  panel.grid(h = 4, v = 4)
  panel.superpose(...)
  },    panel.groups = function(...) {
    panel.histogram(..., border = col)  
    },
  key = list(text = nind.key$text, space = "right", 
    rectangles = list(col = cols[1:4]))
  )

dev.off()

histogram(~dist|method*model.ta, data = qtls[qtls$cd == 1,], 
  groups = nind.t, breaks = 100, 
  col = cols[4:1], xlim = c(0,25), ylim = c(0, 100), 
  xlab = "cM",
  panel = function(...) {
  panel.grid(h = 4, v = 4)
  panel.superpose(...)
  },    panel.groups = function(...) {
    panel.histogram(..., border = col)  
    },
  key = list(text = nind.key$text, space = "right", 
    rectangles = list(col = cols[1:4]))
  )



hs.a = histogram(~dist|nind.t*model.t*method*nmar, data = qtls, breaks = 100)
hs.c = histogram(~dist|nind.t*model.t*method*nmar, data = qtls[qtls$cd == 1,], breaks = 100  )

qtl.a = hs.a$packet.sizes
qtl.c = hs.c$packet.sizes

## copy the packet sizes to excel, sort to make one table.  Use the number of QTL 
## to estimate power and fdr

pc = read.csv("powercalc.csv")
pc.nomd = read.csv("powercalc_nomd.csv")
pc.nomd$method = factor(pc.nomd$method, levels = c("IM", "CIM", "QMS"))
pc.nomd$model = gsub(",", " |", pc.nomd$model)
pc.nomd$model = factor(pc.nomd$model, levels = c("a1 = 5 | a2 = 5",  "a1 = 5 | a2 = 10"))


tiff("Power by nind and model.tiff", width = wi, height = he, 
  units= "px", pointsize = 20)
xyplot(power ~ nind|model, data = pc.nomd, group = method,
    strip = strip.custom(par.strip.text = list(font = 4)), 
    pch = 20, type = "b", lwd = 7, cex = 3 ,
    col = cols[1:3],
    xlab = "Cross Size (N)", ylab = "Power (100 * CD/Total QTL)",
    key = list(text = list(c("IM", "CIM", "QMS"), font = 4),
            lines = list(col = cols[1:3], lwd = 4, pch = 20, cex = 3, type = "b"),
            space = "right")
    )
dev.off()


sapply(c("N = 250", "N = 500", "N = 1000", "N = 1500"), function(IND) 
  sapply(c("IM", "CIM", "QMS"), 
    function(METH) 
    with(qtls[qtls$model == 2 & qtls$nind.t == IND & qtls$method == METH,], 
      cor(residuals(lm(peak ~ detected)), loc)
      )
    )
  )

sapply(c("N = 250", "N = 500", "N = 1000", "N = 1500"), function(IND) 
  sapply(c("IM", "CIM", "QMS"), 
    function(METH) 
    with(qtls[qtls$model == 1 & qtls$nind.t == IND & qtls$method == METH,], 
      cor(peak, loc)
      )
    )
  )

tiff("ID Peak LOD vs distance.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
  
xyplot(peak.lod ~ dist|method*nind.t, groups = model.ta, data = qtls[qtls$cd == 0,], 
  pch = ".", col = cols[c(5,3)],
  xlab = "Distance from Real QTL (cM)", ylab = "Peak LOD Score",
  key = list(text = list(c("a1 =  5 | a2 = 5", "a1 = 5 | a2 = 10"), font = 2, cex = 1.1),
    points = list(col = cols[c(3,5)], pch = 20, cex = 2), space = "right")
  )

dev.off()

tiff("ID Peak LOD vs CI length.tiff", width = wi, height = he, units = "px",
  pointsize = 20)

xyplot(peak.lod ~ ci.length|method*nind.t, groups = model.ta, data = qtls[qtls$cd == 0,], 
  pch = ".", col = cols[c(5, 3)],
  xlab = "Confidence Interval Length (cM)", ylab = "Peak LOD Score",
#auto.key = TRUE
  key = list(text = list(c("a1 =  5 | a2 = 5", "a1 = 5 | a2 = 10"), font = 2, cex = 1.1),
    points = list(col = cols[c(3, 5)], pch = 20, cex = 2), space = "right")
  )
  
dev.off()


par(mfrow = c(3, 3))
sapply(1:9, function(i) hist(rgamma(n = 10000, shape = i, scale = 2), breaks = 1000))
 
tiff("CI length vs Distance by detected QMS.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(ci.length ~ dist|model.t*nind.t, groups = cd, data = qtls[qtls$method == "QMS",], 
  pch = ".", col = cols[1:2],
  cex.axis = 1.6, cex.lab = 2,
  xlab = "Distance from real QTL (cM)", ylab = "Confidence Interval length (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()

tiff("CI length vs Distance by detected CIM.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(ci.length ~ dist|model.t*nind.t, groups = cd, data = qtls[qtls$method == "CIM",], 
  pch = ".", col = cols[1:2],
  cex.axis = 1.6, cex.lab = 2,
  xlab = "Distance from real QTL (cM)", ylab = "Confidence Interval length (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()

tiff("CI length vs Distance by detected IM.tiff", width = wi, height = he, units = "px",
  pointsize = 20)
xyplot(ci.length ~ dist|model.t*nind.t, groups = cd, data = qtls[qtls$method == "IM",], 
  pch = ".", col = cols[1:2],
  cex.axis = 1.6, cex.lab = 2,
  xlab = "Distance from real QTL (cM)", ylab = "Confidence Interval length (cM)",
 key = list(text = list(c("FD", "CD"), font = 2, cex = 1.5),
            points = list(col = cols[1:2], pch = 20, cex = 3), space = "right"
          )
  )
dev.off()


histogram(~dist| method*model.ta, groups = detected, qtls[qtls$cd == 1,], 
  col = cols[c(3,5)],  breaks = 100,
  panel = function(...) {
  panel.grid(h = 4, v = 4)
  panel.superpose(...)
  },    panel.groups = function(...) {
    panel.histogram(..., border = col)  
    }
)
