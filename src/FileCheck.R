##orterun ~/R/bin/R --vanilla <cf2_analyze_q1_world.R > q1a.Rtxt >&1


files.p = list.files(path = "cf2_o", pattern = "RData")
files.d = list.files(path = "cf2_d", pattern = "csv")
## Check first: files.d = files.d[-c(16001,16002)]

files.p = gsub("RData", "csv", files.p)
files.p %in% files.d
files.d %in% files.p

sapply(files.d, function(i) {
 if(i %in% files.p) system(paste("mv ~/CF2/cf2_d/", i, " ~/CF2/cf2_p", sep = ""), intern = TRUE)
 }
 )
