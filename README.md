# feril
Fortran electron repulsion library

@TODO: integrate a faster Boys function from here:
https://aip.scitation.org/doi/suppl/10.1063/5.0062444/suppl_file/dboysfun12.f90.txt
speedups:
	- 1/2 can be \*0.5
	- complex division can be simplified since only `real` part is
	  used


