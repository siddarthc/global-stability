USE_EB = TRUE

ebase = stabilityDriver

CHOMBO_HOME = /work/02138/siddarth/lonestar/CHOMBO/Chombo-3.2/lib
#CHOMBO_HOME = /work/02138/siddarth/CHOMBO/Chombo-3.2/lib
#CHOMBO_HOME = /home/siddarth/CHOMBO/SFD_CHOMBOlib

LibNames = EBAMRTimeDependent EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTimeDependent AMRTools BoxTools
base_dir = .

#compiling both INS and CNS gives compilation errors due to multiple definitions of InflowOutflowIBC
src_dirs = ../libStability/ ../srcStability ../SFD/libSFD ../SFD/srcSFD ../SFD/SFDUtils ../solvers/EBAMRINS/src ../solvers/EBAMRINS/srcLib 
#src_dirs = ../libStability/ ../srcStability ../SFD/libSFD ../SFD/srcSFD ../solvers/EBAMRCNS/src ../solvers/srcSolvers

INPUT = ramp.2d.inputs

include $(CHOMBO_HOME)/mk/Make.example

