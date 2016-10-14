#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "EBExplosionIBCFactory.H"
#include "EBPatchPolytropicFactory.H"
#include "InflowOutflowIBCFactory.H"
#include "EBPatchPolytropic.H"

#include "EBAMRCNSFactory.H"
#include "EBAMRCNS.H"
#include "AMRLevel.H"
#include "AMR.H"
#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "PolyGeom.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "AllRegularService.H"
#include "EBLevelRedist.H"
#include "RedistStencil.H"
#include "SlabService.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBFABView.H"
#include "memtrack.H"
#include "EBAMRCNSParams.H"
#include "GodunovGeom.H"
#include "EBViscousTensorOpFactory.H"
#include "DirichletPoissonEBBC.H"
#include   "NeumannPoissonEBBC.H"
#include "DirichletPoissonDomainBC.H"
#include   "NeumannPoissonDomainBC.H"
#include "DirichletViscousTensorEBBC.H"
#include   "NeumannViscousTensorEBBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include   "NeumannViscousTensorDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include   "NeumannConductivityDomainBC.H"
#include "DirichletConductivityEBBC.H"
#include   "NeumannConductivityEBBC.H"
#include "CH_Attach.H"

#include <iostream>

#include "UsingNamespace.H"

using std::ifstream;
using std::ios;

/***************/
//int InflowOutflowIBC::m_icount = 0; // initialize tO here

void amrGodunov(const Box&      a_domain,
                const RealVect& a_dx)
{
  AMR amr = setupProblem(a_domain, a_dx);
  // run
  int nstop = 0;
  ParmParse ppgodunov;
  ppgodunov.get("max_step",nstop);

  Real stopTime = 0.0;
  ppgodunov.get("max_time",stopTime);

  bool checkForSteadyState = true;
  amr.checkForSteadyState(checkForSteadyState);

  amr.run(stopTime,nstop);

  // output last pltfile and statistics
  //cleanup
  amr.conclude();
//  delete amrg_fact;
}


/***************/
int
main(int a_argc, char* a_argv[])
{

#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
  // registerDebugger();
  // setChomboMPIErrorHandler();
#endif
  {
    // Check for an input file
    char* inFile = NULL;

    EBDebugPoint::s_ivd = IntVect(D_DECL(16,45,0));
    if (a_argc > 1)
      {
        inFile = a_argv[1];
      }
    else
      {
        pout() << "Usage: <executable name> <inputfile>" << endl;
        pout() << "No input file specified" << endl;
        return -1;
      }
    // Parse the command line and the input file (if any)
    ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);
//    ParmParse pp(0, &a_argv[2], NULL, inFile);

    Box coarsestDomain;
    RealVect dx;
    // run amrGodunov
    godunovGeometry(coarsestDomain, dx);

    amrGodunov(coarsestDomain, dx);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  } //scoping trick

#ifdef CH_MPI
  pout() << "dumping timers" << endl;
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
  pout() << "Done." << endl;

}
