#!/bin/bash

for SMRY in S1 MQM ; do 
    for Q in q1 q2 ; do 
	head -n1 ../cf2_o/1000_CF2_${Q}_nind1000_nmar1152_${SMRY}smry.csv | sed -e 's/\"//g' > ../analysis/${SMRY}_${Q}_summary.csv
	for i in 250 500 1000 1500 ; do 
	    cat ../cf2_o/*${Q}_nind${i}*${SMRY}smry.csv | sed -e 's/\"//g' -e '/^lci/d' -e '/^name/d' >> ../analysis/${SMRY}_${Q}_summary.csv
	done
    done
done

