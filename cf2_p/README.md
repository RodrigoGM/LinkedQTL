###cf2_p/

Processed csv cross files.  


After completion of the analysis script, we run system call within in R to move the analyzed cross to this directory.

      ```system(paste("mv ../cf2_d/", cross.i, " ../cf2_p", sep = ""))``` 

A tar ball was created to contain all the files as some file systems crash with this number of files

      ```tar -cvf processed.tar *csv```
