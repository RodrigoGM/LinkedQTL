LinkageAnalysis.cross<-
function(cross, pheno.col=1, filestem='f2 cross', cross.no=cross.no, n.cluster = 2, ...)
{

#cat(paste("Estimating Recombination Fraction of Cross", cross.no, "\n"))
#    cross<- est.rf(cross=cross)
#    plot.rf(cross)

#cat("Calculating Genotype probabilities of Cross", cross.no, "\n")
#    cross<<- calc.genoprob(cross, step=1, error.prob=0.01,
#      map.function="haldane")

  cat("Analyzing cross w/ scanone, cross", cross.no, "\n")
      out.cross<<- scanone(cross, pheno.col=pheno.col,
        method='hk', n.cluster = n.cluster)
        
  cat("...summaryzing scanone, cross", cross.no, "\n")
      smry.out.cross<<-summary.scanone(out.cross, threshold=3.5, format="allpheno")

  cat("Analyzing cross with cim, cross", cross.no, "\n")
      cim.cross<<- mcim(cross, pheno.col=pheno.col, method='hk', window=2,
        error.prob=0.0001, imp.method='imp', n.marcovar=3)

  cat("...summaryzing cim, cross", cross.no, "\n")
      smry.cim.cross<- summary.scanone(cim.cross, threshold=3.5, format="allpheno")

  cat("Estimating LOD confidence intervals from", cross.no, "\n")
    chr<- unique(as.character(smry.out.cross$chr))
     int.p1<- lapply(chr, function(cc) lodint(out.cross, chr=cc, drop=1.5, lodcolumn=1))
     int.p1<-do.call("rbind", int.p1)
     int.p1$cross.no= cross.no
     int.p1$ForPhen=1
     int.p1$nind=250
     int.p1$method="IM"
     
     int.p2<- lapply(chr, function(cc) lodint(out.cross, chr=cc, drop=1.5, lodcolumn=2))
     int.p2<-do.call("rbind", int.p2)
     int.p2$cross.no= cross.no
     int.p2$ForPhen=2
     int.p2$nind=250
     int.p2$method="Scanone"
     
     int.p3<- lapply(chr, function(cc) lodint(out.cross, chr=cc, drop=1.5, lodcolumn=3))
     int.p3<-do.call("rbind", int.p3)
     int.p3$cross.no= cross.no
     int.p3$ForPhen=3
     int.p3$nind=250
     int.p3$method="IM"

     int.p4<- lapply(chr, function(cc) lodint(out.cross, chr=cc, drop=1.5, lodcolumn=4))
     int.p4<-do.call("rbind", int.p4)
     int.p4$cross.no= cross.no
     int.p4$ForPhen=4
     int.p4$nind=250
     int.p4$method="IM"
     
     int.p5<- lapply(chr, function(cc) lodint(out.cross, chr=cc, drop=1.5, lodcolumn=5))
     int.p5<-do.call("rbind", int.p5)
     int.p5$cross.no= cross.no
     int.p5$ForPhen=5
     int.p5$nind=250
     int.p5$method="IM"


lodinterval.s<- rbind(int.p1, int.p2, int.p3, int.p4, int.p5)

    chr<- unique(as.character(smry.cim.cross$chr))
     int.c1<- lapply(chr, function(cc) lodint(cim.cross, chr=cc, drop=1.5, lodcolumn=1))
     int.c1<-do.call("rbind", int.c1)
     int.c1$cross.no= cross.no
     int.c1$ForPhen=1
     int.c1$nind=250
     int.c1$method="CIM"
     
     int.c2<- lapply(chr, function(cc) lodint(cim.cross, chr=cc, drop=1.5, lodcolumn=2))
     int.c2<-do.call("rbind", int.c2)
     int.c2$cross.no= cross.no
     int.c2$ForPhen=2
     int.c2$nind=250
     int.c2$method="CIM"
     
     int.c3<- lapply(chr, function(cc) lodint(cim.cross, chr=cc, drop=1.5, lodcolumn=3))
     int.c3<-do.call("rbind", int.c3)
     int.c3$cross.no= cross.no
     int.c3$ForPhen=3
     int.c3$nind=250
     int.c3$method="CIM"

     int.c4<- lapply(chr, function(cc) lodint(cim.cross, chr=cc, drop=1.5, lodcolumn=4))
     int.c4<-do.call("rbind", int.c4)
     int.c4$cross.no= cross.no
     int.c4$ForPhen=4
     int.c4$nind=250
     int.c4$method="CIM"
     
     int.c5<- lapply(chr, function(cc) lodint(cim.cross, chr=cc, drop=1.5, lodcolumn=5))
     int.c5<-do.call("rbind", int.c5)
     int.c5$cross.no= cross.no
     int.c5$ForPhen=5
     int.c5$nind=250
     int.c5$method="CIM"

lodinterval.cim<- rbind(int.c1, int.c2, int.c3, int.c4, int.c5)


    if (nmar(cross)[1] == 20)  {out.cross.20<<- out.cross
                             cim.cross.20<<- cim.cross }

    if (nmar(cross)[1] == 100) {out.cross.100<<- out.cross
                             cim.cross.100<<- cim.cross}

    if (nmar(cross)[1] == 200) {out.cross.200<<- out.cross
                             cim.cross.200<<- cim.cross}

    if (nmar(cross)[1] == 500) {out.cross.500<<- out.cross
                             cim.cross.500<<- cim.cross}


   cat("Writing CSV files of summaries of cross", cross.no, "\n")
    write.csv(smry.out.cross,
     row.names=FALSE, file=paste(cross.no, "_", filestem, "_S1_Smry.csv", sep=""))

    write.csv(smry.cim.cross, row.names=FALSE,
      file=paste(cross.no, "_", filestem, "_CIM_Smry.csv", sep=""))

   cat("Writing CSV files of Lod Score Intervals of cross", cross.no, "\n")
    write.csv(lodinterval.s, file=paste(cross.no, "_", filestem, "_LOD_Interval_S1.csv", sep=""))
    write.csv(lodinterval.cim, file= paste(cross.no, "_", filestem, "_LOD_Interval_CIM.csv", sep=""))

}


mcim<- function(cross, pheno.col=c(1,2), n.marcovar=3, window=2, method='hk', 
        imp.method='imp', error.prob=0.0001, ...) {
require(qtl) || stop("qtl package not available.")

if (length(pheno.col)<= 1) {stop("you only have one phenotype, please use 'cim'")}

if (length(pheno.col)>= 2) {
        
      cim..<- as.data.frame(lapply(pheno.col, function(x) cim(cross, pheno.col=x, n.marcovar=n.marcovar, 
            window=window, method=method, error.prob=error.prob, ...)))
    
      cim..<- cim..[,c('chr','pos', 'lod', paste('lod', seq(length(pheno.col)-1), sep="."))]
      colnames(cim..)<- c('chr', 'pos', pheno.col)
      }
      
class(cim..)<- c('scanone', 'data.frame')
attr(cim.., "method")= method
attr(cim.., "map")= pull.map(cross)

return(cim..)

} ## mcim
