### cf2_d/
This directory contains un-analyzed simulated cross data in csv format.

Naming convention is as follows :

  * 1      : cross id
  * CF2    : Congenic F2
  * q1|q2  : QTL model (q1 = two qtl of equal effect ; q2 = two qtl of different effect)
  * nind<N>: Number of individuals in the cross
  * nmar<N>: Number of markers on the chromosome

E.G.  1_CF2_q1_nind1000_nmar1152.csv  1_CF2_q2_nind250_nmar96.csv 

Each cross contains two QTL with a 15 cM spacing.  The equal QTL effect model has two additive QTL with effects set as a1 = 7.5, a2 = 7.5 ; whereas the unequal QTL effect model has the qtl with additive effects set as a1 = 5 and a2 = 10.  Dominace effects for all models was 0.
