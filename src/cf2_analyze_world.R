## Analyze 1000*4*4 crosses from 1 chr with 2 QTL spaced at 15, 30 and 45 cM appart
## with varying a and d effect in a semi factorial design
## Setting up a working directory for the QTL effects 1
setwd("./cf2_o")

## error options
options(error = quote(dump.frames("q1_A3", to.file = TRUE)))

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

## maps = load("../C-map.rda")
ls()

system.time({

    np <- 48 ## mpi.universe.size() ##
    cl <- makeCluster(np, type = "SOCK")  ##  type = "MPI")  ##  
    clusterSetupRNG(cl, seed = c(20100905, 19450728, 19750903, 20090830, 19811130, 19811122))
    clusterExport(cl, list = ls())
    clusterEvalQ(cl, {
        require(qtl)
        ls()
    })
    
    
### FILE LISTS ###
    
    cross.file.list <- list.files(path = "../cf2_d", pattern = "CF2_.*csv")
    
    
### PARALLELIZATION ###
    
    parLapply(cl, cross.file.list, function(cross.i) {
        
        ##lapply(cross.file.list, function(cross.i) {
        cat(cross.i, "\n")
        cross.no <- strsplit(cross.i, split = "_")[[1]][1]
        random.sample = readLines("RandomSample.txt")
        ##==============================================================================
        ## 
        cross.a <- read.cross(format = "csv", dir = "../cf2_d", file = cross.i)
        
        cross.a <- est.rf(cross.a)
        cross.a <- calc.genoprob(cross.a)
        phenos <- names(cross.a$pheno)
        
        im.hk.a <- scanone(cross.a, pheno.col = phenos, method = "hk")
        cim.hk.a <- mcim(cross.a, pheno.col = phenos, method = "hk")
        stq.hk.a <- stepwiseqtl.mt(cross.a, pheno.list = phenos, method = "hk", max.qtl = 4)
        
        smry.im.hk.a <- sumarize.scanone(im.hk.a, cross.no = cross.no, cross = cross.a, pheno.list = phenos,
                                         k.dist = 15, mapped.with = "IM")
        
        smry.cim.hk.a <- sumarize.scanone(cim.hk.a, cross.no = cross.no, cross = cross.a, pheno.list = phenos,
                                          k.dist = 15, mapped.with = "CIM")
        
        smry.stq.hk.a <- summary.mtqtl(stq.hk.a, cross.no = cross.no, cross = cross.a, pheno.list = phenos,
                                       k.dist = 15, mapped.with = "stq")
        
        write.csv(rbind(smry.im.hk.a, smry.cim.hk.a),
                  file = paste(gsub("\\.csv", "", cross.i), "_S1smry.csv", sep = ""),
                  row.names = FALSE, append = TRUE)
        
        write.csv(smry.stq.hk.a, file = paste(gsub("\\.csv", "", cross.i), "_MQMsmry.csv", sep = ""),
                  row.names = FALSE, append = TRUE)
        
        save.image(gsub("csv", "rda", cross.i))
        
        if(cross.no %in% random.sample) save(list = ls(pattern = "a"),
                                             file = paste("../cf2_e/", gsub("csv", "RData", cross.i), sep = ""))
        
        system(paste("mv ../cf2_d/", cross.i, " ../cf2_p", sep = ""))
        
        ##==============================================================================
    })
    
})

stopCluster(cl)
