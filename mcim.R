mcim<- function(cross, pheno.col=c(1,2), n.marcovar=3, window=2, method='hk', 
        imp.method='imp', error.prob=0.0001, ...) {
require(qtl) || stop("qtl package not available.")

if (length(pheno.col)<= 1) {stop("you only have one phenotype, please use 'cim'")}

if (length(pheno.col)>= 2) {
        
      out.cim<- as.data.frame(lapply(pheno.col, function(x) cim(cross, pheno.col=x, n.marcovar=n.marcovar, 
            window=window, method=method, error.prob=error.prob)))
    
      out.cim<- out.cim[,c('chr','pos', 'lod', paste('lod', seq(length(pheno.col)-1), sep="."))]
      colnames(out.cim)<- c('chr', 'pos', pheno.col)
      }
      
class(out.cim)<- c('scanone', 'data.frame')
attr(out.cim, "method")= method
attr(out.cim, "map")= pull.map(cross)

return(out.cim)

} ## mcim 


sumarize.scanone<- function(object, cross.no, cross, pheno.list, k.dist, mapped.with){
dd<- data.frame(do.call("rbind",
                        lapply(pheno.list, function(P) {
                          ints<- lodint(object, lodcolumn = grep(P, pheno.list))[ ,c("chr", "pos", P)]
                          ints<- t(ints)
                          colnames(ints)<- c("lci", "peak", "uci")
                          pos<- ints["pos",]
                          lod<- ints[P,]
                          out<- c(pos, lod, P)
                          return(out)
                     })
       ))
colnames(dd) <- c("lci", "peak", "uci", "lci.lod", "peak.lod", "uci.lod", "phenotype")
dd$name <- paste("1@", dd$peak, sep = "")
dd$chr <- 1
dd$peak.marker = find.marker(cross, chr = 1, pos = as.numeric(dd$peak))
dd$cross.no <- cross.no
dd$nmar <- nmar(cross)
dd$nind <- nind(cross)
dd$known.dist <- k.dist
dd$method <- mapped.with


return(dd)
} ## end sumarize.scanone

