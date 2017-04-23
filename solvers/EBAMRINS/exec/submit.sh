#!/bin/bash
#$ -cwd
#$ -l qty.eq.4,xeon
#$ -N CJ_SFD

source /share/cfdlabsetup

#/share/bin/parampich ./IncNavierStokesSolver counter2D_new.xml conditions_counter2D_adaptiveP.xml
#/share/bin/parampich ./IncNavierStokesSolver counter2D_adaptiveSFD.xml.gz counter2D_adaptiveSFD.xml
#/share/bin/parampich ./navierDriver2d.Linux.64.mpicxx.mpifort.MPI.ex sphere.inputs
/share/bin/parampich ./viscousDriver2d.Linux.64.mpicxx.mpif90.MPI.ex sphere.inputs
