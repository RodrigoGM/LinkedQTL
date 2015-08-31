## Generate 1000 crosses from 1 chr with 2 QTL spaced at 15, 30 and 45 cM appart
##  with varying a and d effect in a semi factorial design
##  Setting up a working directory for the QTL effects 1
setwd("./cf2_d")

## error options
options(error = quote(dump.frames("q1", to.file = TRUE)))

.libPaths="/home/ulg/genan/rgularte/lib64/R/library/"
library(qtl, lib.loc = "/home/ulg/genan/rgularte/lib64/R/library/")
library(RColorBrewer, lib.loc = "/home/ulg/genan/rgularte/lib64/R/library/")
library(snow, lib.loc = "/home/ulg/genan/rgularte/lib64/R/library/")
library(rlecuyer, lib.loc = "/home/ulg/genan/rgularte/lib64/R/library/")

sapply(file.path("../src/", c("FileCheck.R", "mcim.R", "SelectMarkers.R", 
                              "sim.cross.f2.npheno.R", "sim.cross.npheno.R", "plotmedians.R",
                              "plotmedians2.R", "LinkageAnalysis.R", "StepwiseAnalysis.R")),
       function(f) {
	 source(paste(f, sep = ""))
       })


load("../data/C-map.rda")
load("../data/seed1.rda")
qtl <- read.csv("../data/CF2-QTL_effects.csv")
## map = sim.map(100, n.mar = 384*4, eq.spacing = TRUE, include.x = FALSE)
## map96 = sim.map(100, n.mar = 96, eq.spacing = TRUE, include.x = FALSE)
## map384 = sim.cross(100, n.mar = 384, eq.spacing = TRUE, include.x = FALSE)
## map1152 = sim.cross(100, n.mar = 1152, eq.spacing = TRUE, include.x = FALSE)
## save(map, map96, map384, map1152, file = "C-map.RData")

set.seed(seed1)

reps <- 1800

Q1 <- data.frame(loc1 = runif(reps, 5, 95-15))
Q1$loc2 <- Q1$loc1 + 15
Q1 <- cbind(Q1, qtl[1,])
write.csv(Q1, file = "CF2-Q1_Model.csv")

q1 <- lapply(1:reps, function(i) rbind(c(Q1[i,"chr"], Q1[i,"loc1"], Q1[i,"a1"], Q1[i,"d1"]),
                                       c(Q1[i,"chr"], Q1[i,"loc2"], Q1[i,"a2"], Q1[i,"d2"])
                                       ))

Q2 <- data.frame(loc1 = runif(reps, 5, 95-15))
Q2$loc2 <- Q2$loc1 + 15
Q2 <- cbind(Q2, qtl[2,])
write.csv(Q2, file = "CF2-Q2_Model.csv")

q2 <- lapply(1:reps, function(i) rbind(c(Q2[i,"chr"], Q2[i,"loc1"], Q2[i,"a1"], Q2[i,"d1"]),
                                       c(Q1[i,"chr"], Q1[i,"loc2"], Q1[i,"a2"], Q1[i,"d2"])
                                       ))


## Function to make all N crosses the actual running is done below
make.world <- function(cl, N = 10, map, nind = 100, qtl.models, qtl., nMar = as.numeric(nmar(map))) {
  parLapply(cl, 1:N, function (i) {
    ##sapply(1:N, function (i) {
      complete.cross <- sim.cross.npheno(map, n.ind = nind, model = qtl.models[[i]],
					 type = "f2p5", map.function = "haldane", error.prob = 0.009, missing.prob = 0.009)
      write.cross(complete.cross, format = "csv", paste(i, "_CF2", "_q", qtl. , "_nind", nind, "_nmar", nMar, sep = ""))
    })    ## end parLapply
    }     ## make.world

Rprof("q1.prof")
system.time({
  
  np <- 48 ## mpi.universe.size() ##
  n.crosses = reps
  ##nMar <-as.numeric(nmar(map))
  maps = list(map, map1152, map384, map96)
  cl <- makeCluster(np, type = "SOCK")  ##  type = "MPI")  ##  
  clusterSetupRNG(cl, seed = c(20100905, 19450728, 19850327, 20090830, 19811130, 19811122))
  clusterExport(cl, list = ls())
  clusterEvalQ(cl, {
      require(qtl)
      ls()
    })
	    
    sapply(maps, function(MAP) {

      make.world(cl = cl, N = n.crosses, map = MAP, nind = 1500, qtl.models = q1, qtl. = 1, nMar = as.numeric(nmar(MAP)))
      make.world(cl = cl, N = n.crosses, map = MAP, nind = 1000, qtl.models = q1, qtl. = 1, nMar = as.numeric(nmar(MAP)))
      make.world(cl = cl, N = n.crosses, map = MAP, nind = 500, qtl.models = q1, qtl. = 1, nMar = as.numeric(nmar(MAP)))
      make.world(cl = cl, N = n.crosses, map = MAP, nind = 250, qtl.models = q1, qtl. = 1, nMar = as.numeric(nmar(MAP)))

      make.world(cl = cl, N = n.crosses, map = MAP, nind = 1500, qtl.models = q2, qtl. = 2, nMar = as.numeric(nmar(MAP)))
      make.world(cl = cl, N = n.crosses, map = MAP, nind = 1000, qtl.models = q2, qtl. = 2, nMar = as.numeric(nmar(MAP)))
      make.world(cl = cl, N = n.crosses, map = MAP, nind = 500, qtl.models = q2, qtl. = 2, nMar = as.numeric(nmar(MAP)))
      make.world(cl = cl, N = n.crosses, map = MAP, nind = 250, qtl.models = q2, qtl. = 2, nMar = as.numeric(nmar(MAP)))
    })
    
    stopCluster(cl)
    
  })
    
    Rprof(NULL)
    
    
