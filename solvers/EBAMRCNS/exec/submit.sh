#!/bin/bash
#$ -cwd
#$ -l qty.eq.8,xeon
#$ -N MPITESTJOB

#/share/bin/parampich ./IncNavierStokesSolver counter2D_new.xml conditions_counter2D_adaptiveP.xml
#/share/bin/parampich ./IncNavierStokesSolver counter2D_adaptiveSFD.xml.gz counter2D_adaptiveSFD.xml
/share/bin/parampich ./navierDriver2d.Linux.64.mpicxx.mpifort.MPI.ex sphere.inputs
#/share/bin/parampich ./navierDriver2d.Linux.64.mpicxx.mpifort.MPI.ex ramp_test.inputs
