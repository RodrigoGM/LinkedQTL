###cf2_p/

Processed csv cross files.  


After completion of the analysis script, we run system call within in R to move the analyzed cross to this directory.

      ```system(paste("mv ../cf2_d/", cross.i, " ../cf2_p", sep = ""))``` 
