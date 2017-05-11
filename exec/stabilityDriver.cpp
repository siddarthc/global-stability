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
#include "EBAMRLinINS.H"
#include "InflowOutflowIBC.H"
#include "InflowOutflowParams.H"
#include "PoiseuilleInflowBCValue.H"
#include "counterJetIBC.H"
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
void setupInflowOutflowIBC(RefCountedPtr<EBIBCFactory>& a_ibc, bool a_homogeneousBC)
{

  // read inputs
  ParmParse pp;

  int flowDir;
  pp.get("flow_dir", flowDir);
  Vector<int> nCells;
  pp.getarr("n_cell",  nCells,0,SpaceDim);

  Real inflowVel = 0.;
  if (!a_homogeneousBC) pp.get("inflow_vel", inflowVel);

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

  bool doAdjoint;
  pp.get("do_adjoint", doAdjoint);

  bool doTransientGrowth;
  pp.get("do_transient_growth", doTransientGrowth);

  bool homogeneousBC = (doAdjoint || doTransientGrowth);

  int orderEBBC = 1;
  pp.query("order_ebbc", orderEBBC);

  bool doPoiseInflow = false;
  pp.query("poiseuille_inflow", doPoiseInflow);

  bool initPoiseData = false;
  pp.query("poiseuille_init", initPoiseData);
  if (initPoiseData)
  pout() << "Doing Poiseuille initialization" << endl;

  RefCountedPtr<PoiseuilleInflowBCValue> poiseBCValue;//make this BaseBCValue if also doing constant values

  bool doWomersleyInflow = false;
  pp.query("womersley_inflow", doWomersleyInflow);

  InflowOutflowParams ibc_params;
  ParseInflowOutflowParams(ibc_params);

  a_ibc = RefCountedPtr<EBIBCFactory>(new InflowOutflowIBCFactory(flowDir,
                                                                  inflowVel,
                                                                  orderEBBC,
                                                                  ibc_params,
                                                                  doSlipWallsHi,
                                                                  doSlipWallsLo,
                                                                  homogeneousBC,
                                                                  doPoiseInflow,
                                                                  initPoiseData,
                                                                  poiseBCValue,
                                                                  doWomersleyInflow));
}
/*********/
void setupCounterJetIBC(RefCountedPtr<EBIBCFactory>& a_ibc, bool a_homogeneousBC)
{

  // read inputs
  ParmParse pp;


  int flowDir;
  pp.get("flow_dir", flowDir);
  Vector<int> nCells;
  pp.getarr("n_cell",  nCells,0,SpaceDim);

  Real jet1inflowVel = 0.;
  if (!a_homogeneousBC) pp.get("jet1_inflow_vel", jet1inflowVel);

  int idoJet2;
  pp.get("do_jet2", idoJet2);
  bool doJet2 = idoJet2==1;

  Real jet2inflowVel = 0.;
  if (!a_homogeneousBC) pp.get("jet2_inflow_vel", jet2inflowVel);

  Real viscosity = 0.0;
  pp.get("viscosity", viscosity);

  int idoSlipWalls;
  pp.get("do_slip_walls", idoSlipWalls);
  bool doSlip = idoSlipWalls==1;
  IntVect doSlipWallsLo = idoSlipWalls*IntVect::Unit;
  IntVect doSlipWallsHi = idoSlipWalls*IntVect::Unit;
  Vector<int> slipWallsLo,slipWallsHi;
  if(doSlip)
    {
      pp.getarr("do_slip_walls_hi",slipWallsHi,0,SpaceDim);
      pp.getarr("do_slip_walls_lo",slipWallsLo,0,SpaceDim);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          doSlipWallsHi[idir] = slipWallsHi[idir];
          doSlipWallsLo[idir] = slipWallsLo[idir];
        }
    }

  bool doAdjoint;
  pp.get("do_adjoint", doAdjoint);

  bool doTransientGrowth;
  pp.get("do_transient_growth", doTransientGrowth);

  bool homogeneousBC = (doAdjoint || doTransientGrowth);

  int orderEBBC = 1;
  pp.query("order_ebbc", orderEBBC);

  bool doJet1PoiseInflow = false;
  pp.query("jet1_poiseuille_inflow", doJet1PoiseInflow);

  bool doJet2PoiseInflow = false;
  pp.query("jet2_poiseuille_inflow", doJet2PoiseInflow);

  bool initPoiseData = false;
  pp.query("poiseuille_init", initPoiseData);
  if(initPoiseData)
  pout() << "Doing Poiseuille initialization" << endl;

  RefCountedPtr<PoiseuilleInflowBCValue> jet1PoiseBCValue;//make this BaseBCValue if also doing constant values
  RefCountedPtr<PoiseuilleInflowBCValue> jet2PoiseBCValue;

  RealVect jet1centerPt, tubeAxis;
  Real jet1tubeRadius;
  Vector<Real> jet1centerPtVect, tubeAxisVect;
  pp.get("jet1_poise_profile_radius", jet1tubeRadius);
  pp.getarr("jet1_poise_profile_center_pt",  jet1centerPtVect,0,SpaceDim);
  pp.getarr("poise_profile_axis",       tubeAxisVect,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    jet1centerPt[idir] = jet1centerPtVect[idir];
    tubeAxis[idir] = tubeAxisVect[idir];
  }

  Real maxVelFactor;//= 1.5 for planar geometry, = 2.0 for cylindrical
  pp.get("poise_maxvel_factor", maxVelFactor);
  Real jet1maxVel = maxVelFactor*jet1inflowVel;

  jet1PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue>(new PoiseuilleInflowBCValue(jet1centerPt,tubeAxis,jet1tubeRadius,jet1inflowVel,flowDir,viscosity));

  RealVect jet2centerPt, jet2TubeOrig, jet2TubeEnd;
  Real jet2tubeRadius;
  Vector<Real> jet2centerPtVect, jet2TubeOrigVect, jet2TubeEndVect;
  pp.get("jet2_poise_profile_radius", jet2tubeRadius);
  pp.getarr("jet2_poise_profile_center_pt",  jet2centerPtVect,0,SpaceDim);
  pp.getarr("jet2_tube_origin",  jet2TubeOrigVect,0,SpaceDim);
  pp.getarr("jet2_tube_end",  jet2TubeEndVect,0,SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    jet2centerPt[idir] = jet2centerPtVect[idir];
    jet2TubeOrig[idir] = jet2TubeOrigVect[idir];
    jet2TubeEnd[idir] = jet2TubeEndVect[idir];
  }

  Real jet2maxVel = maxVelFactor*jet2inflowVel;

  jet2PoiseBCValue = RefCountedPtr<PoiseuilleInflowBCValue>(new PoiseuilleInflowBCValue(jet2centerPt,tubeAxis,jet2tubeRadius,jet2inflowVel,flowDir,viscosity));

  a_ibc = RefCountedPtr<EBIBCFactory>(new counterJetIBCFactory(flowDir,
                                                               doJet2,
                                                               jet1inflowVel,
                                                               jet2inflowVel,
                                                               orderEBBC,
                                                               jet2TubeOrig,
                                                               jet2TubeEnd,
                                                               doSlipWallsHi,
                                                               doSlipWallsLo,
                                                               homogeneousBC,
                                                               doJet1PoiseInflow,
                                                               doJet2PoiseInflow,
                                                               initPoiseData,
							       jet1PoiseBCValue,
                                                               jet2PoiseBCValue));
}
/*********/
void executeStabilityEvaluator(const AMRParameters& a_params,
                               const ProblemDomain& a_coarsestDomain)
{

  CH_TIMERS("stabilityEvaluator");
  CH_TIMER("computeDomainantModes", t3);

  // read inputs
  ParmParse pp;
  
  Real viscosity = 0.0;
  pp.get("viscosity", viscosity);

  bool incOverlapData; 
  pp.get("include_overlap_data", incOverlapData);

  bool doWeighting;
  pp.get("do_weighted_norm", doWeighting);

  bool plotSnapshots;
  pp.get("plot_snapshots", plotSnapshots);

  bool linINS;
  pp.get("do_linearized_INS", linINS);

  bool doAdjoint;
  pp.get("do_adjoint", doAdjoint);

  bool doTransientGrowth;
  pp.get("do_transient_growth", doTransientGrowth);

  bool firstOrderFreDeriv;
  pp.get("do_first_order_fre_deriv", firstOrderFreDeriv);

// begin Baseflow IBC:

  RefCountedPtr<EBIBCFactory> baseflowIBCFact;
//  setupInflowOutflowIBC(baseflowIBCFact, false);
  setupCounterJetIBC(baseflowIBCFact, false);

// end Baseflow IBC

// begin solver IBC for Linearized solver (either direct or adjoint):

  RefCountedPtr<EBIBCFactory> solverIBCFact;

  if (linINS)
  {
//  setupInflowOutflowIBC(solverIBCFact, true);
    setupCounterJetIBC(solverIBCFact, true);
  }

// end solver IBC 

  // get all params
  // get all the parameters for stability solve:
  Real eps, integrationTime;
  std::string baseflowFile;
  pp.get("perturbation_size", eps);
  pp.get("integration_time", integrationTime);
  pp.get("baseflow_file", baseflowFile);
  
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

  // initial data
  bool initWithRandomData;
  pp.get("init_random_data", initWithRandomData);
  Vector<string> initialDataFile;
  if (!initWithRandomData)
  {
    pp.getarr("initial_data_file",  initialDataFile,0,numBlocks-1);
  }

  RefCountedPtr<EBAMRINSInterfaceFactory> INSFact = RefCountedPtr<EBAMRINSInterfaceFactory>(new EBAMRINSInterfaceFactory(a_params, baseflowIBCFact, solverIBCFact, a_coarsestDomain, viscosity, plotSnapshots, linINS, doAdjoint, doTransientGrowth, firstOrderFreDeriv));

  RefCountedPtr<ChomboSolverInterfaceFactory> solverFact = static_cast<RefCountedPtr<ChomboSolverInterfaceFactory> >(INSFact);

  Teuchos::RCP<TrilinosChomboInterfaceFactory> ChomboFact = Teuchos::rcp (new TrilinosChomboInterfaceFactory(solverFact, initialDataFile, incOverlapData, doWeighting));

  Teuchos::RCP<TrilinosSolverInterfaceFactory> castChomboFact = static_cast<Teuchos::RCP<TrilinosSolverInterfaceFactory> >(ChomboFact);

  MPI_Comm* wcomm = NULL;

#ifdef CH_MPI
  wcomm = &(Chombo_MPI::comm);
#endif

  StabilityEvaluator stabEval(eps, integrationTime, baseflowFile, castChomboFact, wcomm);

  CH_START(t3);
  stabEval.computeDominantModes(evTol, nev, numBlocks, blockSize, maxRestarts, sortEV, false, true, plotEVComps, initWithRandomData);
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
