SelectMarkers2Drop<- 
function(cross, map, n.markers.chr = 200, keep.markers.chr = 20, ...) 
{

bins<<- cbind(seq(from = 1, to = n.markers.chr, by = floor(n.markers.chr/keep.markers.chr)),
         seq(from = ceiling(n.markers.chr/keep.markers.chr), to = n.markers.chr, 
            by = (n.markers.chr/keep.markers.chr)))
  (dim(bins))

marker.bins <- sapply(map, function(chr) apply(bins, 1, function(bin) list(names(chr)[bin[1]:bin[2]])))
                                                              
all.markers = as.character(unlist(lapply(map, names)))

markers.to.drop<- all.markers[all.markers %in% as.character(unlist(apply(marker.bins, 2, function(x) lapply(x, function(y) lapply(y, function(i) sample(i,1)))))) == FALSE]

cross<- drop.markers(cross, markers.to.drop)

cross
}
