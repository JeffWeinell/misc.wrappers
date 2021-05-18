//Number of population samples (demes)
3
//Population effective sizes (number of genes)
N0$
N1$
N2$
//Sample sizes
7
7
10
//Growth rates: negative growth implies population expansion
0
0
0
//Number of migration matrices : 0 (or 1 and the zero matrix) implies no migration between demes
2
//migrationmatrix
0 0 0
0 0 0
0 0 0
//migrationmatrix
0 MIG1$ MIG1$
MIG1$ 0 MIG1$
MIG1$ MIG1$ 0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
2 historical event
Tmig1$ 0 0 0 1 0 1
Tdiv1$ 0 1 1 1 0 0
Tdiv2$ 1 2 1 1 0 0
//Number of independent loci [chromosome]
4176
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loc., rec. rate, and mut rate + optional paramters
SNP 1 0 0 0.01
