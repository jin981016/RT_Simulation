module purge
module load intel/19.1.3 impi/2019.9
mpiifort -g -traceback -heap-arrays 200 -check bounds -ipo -O3 -no-prec-div -fp-model fast=2 -cpp -DMPI -c *_mod.f90
mpiifort -g -traceback -heap-arrays 200 -check bounds -ipo -O3 -no-prec-div -fp-model fast=2 -cpp -DMPI -o CIV_12_14_1_50.out main_12_14_1_50.f90 *.o



#mpiifort -cpp -DMPI -c *_mod.f90
#mpiifort -cpp -DMPI -o CIV.out main_outflow.f90 *.o
