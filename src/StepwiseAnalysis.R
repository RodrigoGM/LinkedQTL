#require(snow)||stop("package 'snow' is not available")
#cl<- makeSOCKcluster(n.cluster)
#clusterStopped<- FALSE
#on.exit(if(!clusterStopped) stopCluster(cl))
#clusterSetupRNG(cl)
#clusterEvalQ(cl, require(qtl))
## FUN()
#clusterStopped<- TRUE
#stopCluster(cl)

stepwiseqtl.mt<- function(cross, pheno.list=1:5, cross.no, method, max.qtl) {

if(missing(pheno.list)) stop("you need to provid the names of the phenotypes of interest")

out<-  lapply(pheno.list, function(PHEN) {
      cat("\n\n\nStepwiseqtl for phenotype", PHEN, "\n\n\n")
      stepwiseqtl(cross, pheno.col = PHEN, method = "hk", max.qtl = max.qtl,
      additive.only = TRUE, keeplodprofile = TRUE, keeptrace = TRUE)})
      names(out)<- pheno.list
class(out) <- c("mtqtl", "list")
return(out)

}  ## end stepwiseqtl.mt


summary.mtqtl<- function(mtqtl, cross.no, cross, pheno.list, k.dist, mapped.with) {
  cat("...summaryzing stq, cross", cross.no, "\n")
  
  if(missing(pheno.list)) pheno.list = names(mtqtl)
  
qtlFound<- do.call("rbind", 
  lapply(pheno.list, function(PHENO) {
    tmp <- as.data.frame(summary(mtqtl[[PHENO]]))
    tmp$phenotype <- as.character(PHENO)
    tmp$pLOD <- as.character(attr(mtqtl[[PHENO]], "pLOD"))
    tmp$n.qtl <- as.character(mtqtl[[PHENO]]$n.qtl)
  tmp
  })
)  # qtlFound

qtlFound$peak.marker<- find.marker(cross, chr = qtlFound$chr, pos = qtlFound$pos)

ci.stq <- do.call("rbind", lapply(mtqtl, function(STQ) 
  do.call("rbind", lapply(1:STQ$n.qtl, function(QI) lodint(STQ, qtl.index = QI)$pos))
  )) # end CI placment
  
rownames(ci.stq) = rownames(qtlFound)
colnames(ci.stq) = c("lci", "peak", "uci")

lodProfiles<- lapply(pheno.list, function(PHEN) {
    lodprofile = attr(mtqtl[[PHEN]], "lodprofile")
    names(lodprofile) = NULL
    lodprofile = do.call("rbind", lodprofile)
}) # end LOD placement
names(lodProfiles) <- pheno.list

qtlFound$lod <- apply(qtlFound, 1, function(i) {
   lod.tmp = lodProfiles[[i[5]]]
   lod.tmp[i["peak.marker"], "lod"]
})

qtlFound <- cbind(qtlFound, ci.stq)

qtlFound$cross.no <- as.numeric(cross.no)
qtlFound$nmar <- as.numeric(nmar(cross))[1]
qtlFound$nind <- as.numeric(nind(cross))
qtlFound$k.dist <- as.numeric(k.dist)
qtlFound$method <- mapped.with

cat("Stepwiseqtl Analysis for Cross", cross.no, "is complete\n \n")

class(qtlFound) <- c("data.frame")
return(qtlFound)
}## end summary.mtqtl
