gcc -c GA.c
 
gcc -c  read.c
 
g++   -I/home/ashwin/lammps-6Apr13/src/STUBS -I/home/ashwin/lammps-6Apr13/src -L/home/ashwin/lammps-6Apr13/src /home/ashwin/lammps-6Apr13/src/STUBS/libmpi_stubs.a  -llammps_serial -lfftw3 -lm -lstdc++   -lpthread -c lammps.cpp #

gcc read.o  Lammps.o GA.o -lstdc++   -I/home/ashwin/lammps-6Apr13/src/STUBS -I/home/ashwin/lammps-6Apr13/src -L/home/ashwin/lammps-6Apr13/src /home/ashwin/lammps-6Apr13/src/STUBS/libmpi_stubs.a  -llammps_serial -lfftw3 -lm -lstdc++  -lrt  -lpthread-#

