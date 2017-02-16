/*
 *
 *
 *
 *
 */

#include "ParmParse.H"
#include "CH_HDF5.H"
#include "parstream.H"

#include "BaseIVFactory.H"
#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "ParmParse.H"
#include "BRMeshRefine.H"
#include "LoadBalance.H"
#include "GeometryShop.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBCellFactory.H"
#include "VoFIterator.H"
#include "EBArith.H"
#include "DebugDump.H"
#include "EBDebugDump.H"
#include "EBFABView.H"
#include "memtrack.H"
#include "AMRINSUtils.H"
#include "CH_Attach.H"
#include "EBAMRNoSubcycle.H"
#include "InflowOutflowIBC.H"
#include "InflowOutflowParams.H"
#include "PoiseuilleInflowBCValue.H"
#include "EBFABView.H"
#include <iostream>
#include "memusage.H"
#include "memtrack.H"

#include "UsingNamespace.H"

#include "StabilityEvaluator.H"
#include "TrilinosChomboInterfaceFactory.H"
#include "EBAMRINSInterfaceFactory.H"

// Trilinos library includes
#include "TrilinosSolverInterfaceFactory.H"
#include <Teuchos_RCPDecl.hpp>
#ifdef CH_MPI
#include <Epetra_MpiComm.h>
#else
#include <Eperta_SerialComm.h>
#endif

/*********/
void executeStabilityEvaluator(const AMRParameters& a_params,
                               const ProblemDomain& a_coarsestDomain)
{

  CH_TIMERS("stabilityEvaluator");
  CH_TIMER("computeDomainantModes", t3);

  // read inputs
  ParmParse pp;
  
  int flowDir;
  pp.get("flow_dir", flowDir);
  Vector<int> nCells;
  pp.getarr("n_cell",  nCells,0,SpaceDim);

  // hardwiring the inflow velocity to 0 because it has to be 0 for stability calculations
  Real inflowVel;
  pp.get("inflow_vel", inflowVel);

  Real viscosity = 0.0;
  pp.get("viscosity", viscosity);

  int idoSlipWalls;
  pp.get("do_slip_walls", idoSlipWalls);
  bool doSlip = idoSlipWalls==1;
  IntVect doSlipWallsLo = idoSlipWalls*IntVect::Unit;
  IntVect doSlipWallsHi = idoSlipWalls*IntVect::Unit;
  Vector<int> slipWallsLo,slipWallsHi;
  if (doSlip)
    {
      pp.getarr("do_slip_walls_hi",slipWallsHi,0,SpaceDim);
      pp.getarr("do_slip_walls_lo",slipWallsLo,0,SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          doSlipWallsHi[idir] = slipWallsHi[idir];
          doSlipWallsLo[idir] = slipWallsLo[idir];
        }
    }

  int orderEBBC = 1;
  pp.query("order_ebbc", orderEBBC); 

  InflowOutflowParams ibc_params;
  ParseInflowOutflowParams(ibc_params);

  bool incOverlapData; 
  pp.get("include_overlap_data", incOverlapData);

  bool doWeighting;
  pp.get("do_weighted_norm", doWeighting);

  bool plotSnapshots;
  pp.get("plot_snapshots", plotSnapshots);

  bool linINS;
  pp.get("do_linearized_INS", linINS);

  bool firstOrderFreDeriv;
  pp.get("do_first_order_fre_deriv", firstOrderFreDeriv);

//  InflowOutflowIBCFactory ibc(flowDir, inflowVel, orderEBBC, ibc_params, doSlipWallsHi, doSlipWallsLo);

  RefCountedPtr<InflowOutflowIBCFactory> ibcFact= RefCountedPtr<InflowOutflowIBCFactory>(new InflowOutflowIBCFactory(flowDir, inflowVel, orderEBBC, ibc_params, doSlipWallsHi, doSlipWallsLo));

  RefCountedPtr<EBIBCFactory> castIBCFact = static_cast<RefCountedPtr<EBIBCFactory> >(ibcFact);

  RefCountedPtr<EBAMRINSInterfaceFactory> INSFact = RefCountedPtr<EBAMRINSInterfaceFactory>(new EBAMRINSInterfaceFactory(a_params, castIBCFact, a_coarsestDomain, viscosity, plotSnapshots, linINS, firstOrderFreDeriv));

  RefCountedPtr<ChomboSolverInterfaceFactory> solverFact = static_cast<RefCountedPtr<ChomboSolverInterfaceFactory> >(INSFact);

  Teuchos::RCP<TrilinosChomboInterfaceFactory> ChomboFact = Teuchos::rcp (new TrilinosChomboInterfaceFactory(solverFact, incOverlapData, doWeighting));

  Teuchos::RCP<TrilinosSolverInterfaceFactory> castChomboFact = static_cast<Teuchos::RCP<TrilinosSolverInterfaceFactory> >(ChomboFact);

  // get all the parameters for stability solve:
  Real eps, integrationTime;
  std::string baseflowFile;
  pp.get("perturbation_size", eps); 
  pp.get("integration_time", integrationTime);
  pp.get("baseflow_file", baseflowFile);
  
  MPI_Comm* wcomm = NULL;

#ifdef CH_MPI
  wcomm = &(Chombo_MPI::comm);
#endif

  StabilityEvaluator stabEval(eps, integrationTime, baseflowFile, castChomboFact, wcomm);

  // get parameters for computing dominant modes
  Real evTol;
  int nev, numBlocks, blockSize, maxRestarts;
  std::string sortEV;
  pp.get("eigenvalue_tol", evTol);
  pp.get("num_eigenValues", nev);
  pp.get("num_blocks", numBlocks);
  pp.get("block_size", blockSize);
  pp.get("max_restarts", maxRestarts);
  pp.get("sort_EV", sortEV);

  int nplotEVComps;
  pp.get("plot_num_EVComps", nplotEVComps);
  std::vector<int> plotEVComps;
  pp.getarr("plot_EVComps",plotEVComps,0,nplotEVComps);
  

  CH_START(t3);
  stabEval.computeDominantModes(evTol, nev, numBlocks, blockSize, maxRestarts, sortEV, false, true, plotEVComps);
   
}

/*********/
int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif

  { //scoping trick

    CH_TIMERS("stabilityRun_timers");
    CH_TIMER("define_geometry", t1);
    CH_TIMER("run", t2);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    // Check for an input file
    char* inFile = NULL;

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

    ProblemDomain coarsestDomain;
    RealVect coarsestDx;

    // specializing for INS:
    AMRParameters params;
    getAMRINSParameters(params, coarsestDomain);

    int numFilt;
    pp.get("num_filter_iterations", numFilt);
    params.m_numFilterIterations = numFilt;

    int gphiIterations;
    pp.get("num_gphi_iterations", gphiIterations);
    params.m_gphiIterations = gphiIterations;

    int initIterations;
    pp.get("num_init_iterations", initIterations);
    params.m_initIterations = initIterations;

    bool doRegridSmoothing;
    pp.get("do_regrid_smoothing", doRegridSmoothing);
    params.m_doRegridSmoothing = doRegridSmoothing;

    CH_START(t1);
    //define geometry
    AMRINSGeometry(params, coarsestDomain);
    CH_STOP(t1); 

    CH_START(t2);
    executeStabilityEvaluator(params, coarsestDomain);
    CH_STOP(t2);

    EBIndexSpace* ebisPtr = Chombo_EBIS::instance();
    ebisPtr->clear();

  } // end scoping

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

}
