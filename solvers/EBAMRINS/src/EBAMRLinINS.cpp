/*
 *
 *
 *
 *
 */

#include "EBAMRLinINS.H"
#include "ParmParse.H"
#include "EBAMRNoSubcycleF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "MeshRefine.H"
#include "BRMeshRefine.H"
#include "EBEllipticLoadBalance.H"
#include "EBArith.H"
#include "EBPWLFineInterp.H"
#include "EBCoarseAverage.H"
#include "EBFluxFactory.H"
#include "EBCellFactory.H"
#include "EBLevelAdvect.H"
#include "EBGradDivFilter.H"
#include "EBPatchAdvect.H"
#include "REAL.H"
#include "EBPhysIBCFactory.H"
#include "EBAMRIO.H"
#include "BaseIFFactory.H"
#include "EBLevelRedist.H"
#include "BaseIVFactory.H"
#include "EBConstantCFInterp.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "NeumannPoissonEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "InflowOutflowIBC.H"
#include "EBNormalizeByVolumeFraction.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "memusage.H"
#include "memtrack.H"

extern Real g_simulationTime;
#define debugIV IntVect(D_DECL(16, 3, 0))

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/
EBAMRLinINS::
EBAMRLinINS(const AMRParameters&      a_params,
            const EBIBCFactory&       a_baseflowIBCFact,
            const EBIBCFactory&       a_solverIBCFact,
            const ProblemDomain&      a_coarsestDomain,
            Real                      a_viscosity,
            bool                      a_doAdjoint,
            const EBIndexSpace* const a_ebisPtr):
  m_ebisPtr(a_ebisPtr)
{
  if (a_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::EBAMRLinINS" << endl;
    }

  //set parameters of the run
  m_params    = a_params;
  m_viscosity = a_viscosity;
  m_viscousCalc = (m_viscosity > 0);

  //create initial and boundary condition object
  m_baseflowIBC    = a_baseflowIBCFact.create();
  m_solverIBC      = a_solverIBCFact.create(); 

  m_doAdjoint = a_doAdjoint;

  //resize vectors and set them where we can
  Real coarsestDx = m_params.m_domainLength/Real(a_coarsestDomain.size(0));
  int nlevels = m_params.m_maxLevel + 1;
  m_domain.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);

  m_ebisl.resize(nlevels);
  m_eblg.resize(nlevels);

  m_quadCFI.resize(nlevels);
  m_aveOper.resize(nlevels);
  m_aveSpac.resize(nlevels);
  m_ebLevAd.resize(nlevels);
  m_fluxReg.resize(nlevels);
  m_velo.resize(nlevels, NULL);
  m_pres.resize(nlevels, NULL);
  m_gphi.resize(nlevels, NULL);

  m_coveredFaceLitLo.resize(nlevels, NULL);
  m_coveredFaceLitHi.resize(nlevels, NULL);
  m_coveredSetsLitLo.resize(nlevels, NULL);
  m_coveredSetsLitHi.resize(nlevels, NULL);

  // base flow data
  m_baseVelo.resize(nlevels, NULL);
  m_baseAdvVelo.resize(nlevels, NULL);
  m_coveredBaseAdvVelLo.resize(nlevels, NULL);
  m_coveredBaseAdvVelHi.resize(nlevels, NULL);

  allocateDataHolders();

  m_domain[0] = a_coarsestDomain;
  m_dx[0]     =   coarsestDx;
  for (int ilev = 1; ilev < nlevels; ilev++)
    {
      CH_assert(m_params.m_refRatio[ilev-1] > 0);
      m_domain[ilev] = refine(m_domain[ilev-1], m_params.m_refRatio[ilev-1]);
      m_dx[ilev] = m_dx[ilev-1]/Real(m_params.m_refRatio[ilev-1]);
    }
  m_prescribedDt = -1.0;
  m_useFixedDt = false;
  m_steadyState = false;
  m_stopAdvance = false;

  m_ccProjector  = NULL;
  m_macProjector = NULL;
  m_time = 0.0;
  m_curStep = 0;
  m_dt = -1.0;

  //setup still needs to get called
  m_isSetup  = false;
  m_pointsUpdated = 0;
}
/*********/
void EBAMRLinINS::
allocateDataHolders()
{
  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      m_velo[ilev]   = new LevelData<EBCellFAB>();
      m_pres[ilev]   = new LevelData<EBCellFAB>();
      m_gphi[ilev]   = new LevelData<EBCellFAB>(); 

      m_coveredFaceLitLo[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredFaceLitHi[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredSetsLitLo[ilev] = new LayoutData< Vector< IntVectSet > >      ();
      m_coveredSetsLitHi[ilev] = new LayoutData< Vector< IntVectSet > >      ();

      m_baseVelo[ilev] = new LevelData<EBCellFAB>();
      m_baseAdvVelo[ilev] = new LevelData<EBFluxFAB>();

      m_coveredBaseAdvVelLo[ilev]  = new LayoutData<Vector<BaseIVFAB<Real> * > >();
      m_coveredBaseAdvVelHi[ilev]  = new LayoutData<Vector<BaseIVFAB<Real> * > >();
    }
}
/*********/
EBAMRLinINS::
~EBAMRLinINS()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::~EBAMRLinINS" << endl;
    }
  delete m_baseflowIBC;
  delete m_solverIBC;

  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      delete m_velo[ilev];
      delete m_pres[ilev];
      delete m_gphi[ilev];

      delete m_coveredFaceLitLo[ilev];
      delete m_coveredFaceLitHi[ilev];
      delete m_coveredSetsLitLo[ilev];
      delete m_coveredSetsLitHi[ilev];

      delete m_baseVelo[ilev];
      delete m_baseAdvVelo[ilev];
  
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*m_coveredBaseAdvVelLo[ilev])[dit()][idir];
              delete (*m_coveredBaseAdvVelHi[ilev])[dit()][idir];
            }
        }
      delete m_coveredBaseAdvVelLo[ilev];
      delete m_coveredBaseAdvVelHi[ilev];
    }

  for (int ilev=0; ilev<m_params.m_maxLevel; ilev++)
    {
      m_baseAdvVelo[ilev]          = NULL;
      m_coveredBaseAdvVelLo[ilev]  = NULL;
      m_coveredBaseAdvVelHi[ilev]  = NULL;
    }

  if (m_ccProjector !=  NULL)
    {
      delete m_ccProjector;
    }
  if (m_macProjector !=  NULL)
    {
      delete m_macProjector;
    }
}
/*********/
void EBAMRLinINS::
useFixedDt(Real a_dt)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::useFixedDt" << endl;
    }
  m_dt = a_dt;
  m_prescribedDt = a_dt;
  m_useFixedDt = true;
}
/*********/
void EBAMRLinINS::
defineEBISLs()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::defineEBISLs" << endl;
    }
  int numEBGhost = 6;//number of ghost cells in EBISL
  m_eblg.resize(m_finestLevel+1);
  //now define all of the storage we need
  RefCountedPtr<EBPhysIBCFactory> advectBC = m_solverIBC->getVelAdvectBC(0); //this gets reset
  RefCountedPtr<EBPatchAdvectFactory> fact = RefCountedPtr<EBPatchAdvectFactory> (new EBPatchAdvectFactory(advectBC, m_params.m_useLimiting));

  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], numEBGhost, m_ebisPtr);
      m_ebisl[ilev] = m_eblg[ilev].getEBISL();
      DisjointBoxLayout coarDBL;
      EBISLayout        coarEBISL;
      int refRat = 2;
       if (ilev > 0)
         {
           coarDBL =  m_grids[ilev-1];
           coarEBISL= m_ebisl[ilev-1];
           refRat = m_params.m_refRatio[ilev-1];
         }
       bool hasCoarser = (ilev > 0);
       bool hasFiner   = (ilev < m_finestLevel);
       m_ebLevAd[ilev]  = RefCountedPtr<EBLevelAdvect>(new EBLevelAdvect(m_grids[ilev],
                                                                         coarDBL,
                                                                         m_ebisl[ilev],
                                                                         coarEBISL,
                                                                         ProblemDomain(m_domain[ilev]),
                                                                         refRat,
                                                                         m_dx[ilev]*RealVect::Unit,
                                                                         hasCoarser,
                                                                         hasFiner,
                                                                         &(*fact),
                                                                         m_ebisPtr));

       /**/
       /**/
      if (ilev > 0)
        {
          //always one component for quadcfi---only way to get reuse
          int nvarquad = 1;
          m_quadCFI[ilev]  = RefCountedPtr<EBQuadCFInterp>(new  EBQuadCFInterp(m_grids[ilev], m_grids[ilev-1],
                                                                               m_ebisl[ilev], m_ebisl[ilev-1],
                                                                               m_domain[ilev-1],
                                                                               m_params.m_refRatio[ilev-1], nvarquad,
                                                                               *m_eblg[ilev].getCFIVS(),
                                                                               m_ebisPtr ));

          m_aveOper[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(m_grids[ilev], m_grids[ilev-1],
                                                                                m_ebisl[ilev], m_ebisl[ilev-1],
                                                                                m_domain[ilev-1],
                                                                                m_params.m_refRatio[ilev-1], nvarquad,
                                                                                m_ebisPtr));

          m_aveSpac[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(m_grids[ilev], m_grids[ilev-1],
                                                                                m_ebisl[ilev], m_ebisl[ilev-1],
                                                                                m_domain[ilev-1],
                                                                                m_params.m_refRatio[ilev-1], SpaceDim,
                                                                                m_ebisPtr));
        }
    }
  for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
      if (ilev < m_finestLevel)
        {
          m_fluxReg[ilev]  = RefCountedPtr<EBFastFR>(new EBFastFR(m_eblg[ilev+1], m_eblg[ilev], m_params.m_refRatio[ilev], SpaceDim));
        }
    }
  long long totalPoints = 0;
  long long totalBoxes  = 0;
  int numLevels = m_finestLevel + 1;
  for (int ilev = 0; ilev < numLevels; ilev++)
    {
      long long pointsThisLevel = 0;
      for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
        {
          pointsThisLevel += m_grids[ilev][lit()].numPts();
        }
      totalPoints += pointsThisLevel;
      totalBoxes += m_grids[ilev].size();
      pout() << "getAllIrregRefineLayouts:level[" << ilev
             << "], number of boxes = " << m_grids[ilev].size()
             << ", number of points = " << pointsThisLevel << endl;
    }
  pout() << "getAllIrregRefineLayouts:"
         <<  "   total boxes = " << totalBoxes
         <<  ", total points = " << totalPoints <<  endl;
}
/*********/
void EBAMRLinINS::
defineNewVel(const int a_startLevel)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::defineNewVel" << endl;
    }

  int startLevel = 0;
  if (a_startLevel > startLevel)
    {
      startLevel=a_startLevel;
    }
  for (int ilev = startLevel; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_velo[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
      m_baseVelo[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      m_baseAdvVelo[ilev]->define(m_grids[ilev], 1,  3*IntVect::Unit, ebfluxfact);
    }
  EBAMRDataOps::setToZero(m_velo);
  EBAMRDataOps::setToZero(m_baseVelo);
  EBAMRDataOps::setToZero(m_baseAdvVelo);
}
/*********/
void EBAMRLinINS::
definePressure(const int a_startLevel)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::definePressure" << endl;
    }
  int startLevel = 0;
  if (a_startLevel > startLevel)
    {
      startLevel=a_startLevel;
    }
  for (int ilev = startLevel; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_gphi[ilev]->define(m_grids[ilev], SpaceDim,    IntVect::Zero, ebcellfact);
      m_pres[ilev]->define(m_grids[ilev],        1,    IntVect::Zero, ebcellfact);
    }
  EBAMRDataOps::setToZero(m_gphi);
  EBAMRDataOps::setToZero(m_pres);
}
/*********/
void EBAMRLinINS::
defineProjections()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::defineProjections" << endl;
    }
  if (m_ccProjector != NULL)
    {
      delete m_ccProjector;
      m_ccProjector = NULL;
    }
  if (m_macProjector != NULL)
    {
      delete m_macProjector;
      m_macProjector = NULL;
    }
  RefCountedPtr<BaseDomainBCFactory> macBCVel = m_solverIBC->getMACVelBC();
  RefCountedPtr<BaseDomainBCFactory> celBCPhi = m_solverIBC->getPressBC();
  RefCountedPtr<BaseEBBCFactory>     ebbcVelo = m_solverIBC->getVelocityEBBC(0);
  RefCountedPtr<BaseEBBCFactory>     ebbcPhi  = m_solverIBC->getPressureEBBC();

  int numLevels = m_finestLevel + 1;
  RealVect coarDxVect = m_dx[0]*RealVect::Unit;

  ParmParse pp;

  pp.query("mg_num_smooths", m_params.m_numSmooth);

  pp.query("mg_relax_type", m_params.m_relaxType);
  bool lazy = false;
  pp.query("mg_relax_lazy", lazy);

  int bottomSolverType = 0;
  pp.query("mg_bottom_solver", bottomSolverType);

  int numPreCond = 4;
  pp.query("mg_num_precond", numPreCond);

  pp.query("mg_hang", m_params.m_hang);
  pp.query("mg_tolerance", m_params.m_tolerance);

  EBAMRPoissonOp::doLazyRelax(lazy);

  int mgCoarsenLimit = 2;
  pp.query("mg_coarsen_limit", mgCoarsenLimit);
  EBAMRPoissonOpFactory::setTestRef(mgCoarsenLimit);
  int whichReflux = 0;
  EBAMRPoissonOpFactory::setWhichReflux(whichReflux);

  m_macProjector =   new EBCompositeMACProjector(m_eblg, m_params.m_refRatio, m_quadCFI,
                                                 coarDxVect, RealVect::Zero,
                                                 macBCVel, celBCPhi, ebbcPhi,
                                                 m_params.m_subtractOffMean, numLevels,
                                                 m_params.m_verbosity,numPreCond,m_time,m_params.m_relaxType,bottomSolverType);

  m_macProjector->setSolverParams(m_params.m_numSmooth, m_params.m_iterMax, m_params.m_mgCycle,
                                  m_params.m_hang, m_params.m_tolerance,
                                  m_params.m_verbosity);

  m_ccProjector  =   new EBCompositeCCProjector( m_eblg, m_params.m_refRatio, m_quadCFI,
                                                 coarDxVect, RealVect::Zero,
                                                 macBCVel, celBCPhi, ebbcPhi,
                                                 m_params.m_subtractOffMean, numLevels,
                                                 m_params.m_verbosity, numPreCond, m_time, m_params.m_relaxType,
                                                 bottomSolverType,m_macProjector);
  if (m_viscousCalc)
    {
      if (m_params.m_orderTimeIntegration == 2)
        {
          pout() << "using TGA for viscous solver" << endl;
        }
      else if (m_params.m_orderTimeIntegration == 1)
        {
          pout() << "using backward Euler for viscous solver" << endl;
        }
      else
        {
          MayDay::Error("EBAMRLinINS::defineProjections -- bad order time integration");
        }

      int lbase = 0;
      int lmax = m_finestLevel;

      LinearSolver<LevelData<EBCellFAB> >* bottomSolverPtr = NULL;
      if (bottomSolverType == 0)//BiCGStab
        {
          bottomSolverPtr = &m_bottomSolver;
        }
      else if (bottomSolverType == 1)//simple relaxation
        {
          bottomSolverPtr = &m_bottomSolverSimp;
        }
      else
        {
          MayDay::Error("EBAMRLinINS::defineProjections (1) -- bad bottom solver type");
        }

      Vector<LevelData<EBCellFAB>* > phi(numLevels);
      Vector<LevelData<EBCellFAB>* > rhs(numLevels);
      for (int ilev = 0; ilev < numLevels; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          phi[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1,   3*IntVect::Unit, ebcellfact);
          rhs[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1,   3*IntVect::Unit, ebcellfact);
        }

      m_bottomSolver.m_verbosity = m_params.m_verbosity;
      int numBotSmooth = 64;
      pp.query("mg_num_bottom_smooths", numBotSmooth);
      m_params.m_numBotSmooth = numBotSmooth;
      m_bottomSolverSimp.setNumSmooths(m_params.m_numBotSmooth*m_params.m_numSmooth);
      pp.query("mg_num_cycles", m_params.m_mgCycle);
      pp.query("mg_iter_max", m_params.m_iterMax);
      pp.query("mg_norm_thresh", m_params.m_normThresh);
      int numLevels = m_finestLevel + 1;

      Real alpha = 1.0;
      Real beta = m_viscosity;

      ProblemDomain level0Dom = m_eblg[0].getDomain();

      for (int idir = 0;  idir < SpaceDim; idir++)
        {
          DirichletPoissonEBBC::s_velComp = idir;
          RefCountedPtr<BaseDomainBCFactory> celBCVel = m_solverIBC->getVelBC(idir);

          EBAMRPoissonOpFactory opFactory(m_eblg,
                                          m_params.m_refRatio,
                                          m_quadCFI,
                                          coarDxVect,
                                          RealVect::Zero,
                                          numPreCond,
                                          m_params.m_relaxType,
                                          celBCVel,
                                          ebbcVelo,
                                          alpha,
                                          beta,
                                          m_time,
                                          3*IntVect::Unit,
                                          3*IntVect::Unit);

          m_solver[idir] = RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > > (new AMRMultiGrid<LevelData<EBCellFAB> >);
          m_solver[idir]->define(level0Dom,
                                 opFactory,
                                 bottomSolverPtr,
                                 numLevels);

          m_solver[idir]->m_pre        =  m_params.m_numSmooth;
          m_solver[idir]->m_post       =  m_params.m_numSmooth;
          m_solver[idir]->m_bottom     =  m_params.m_numSmooth;
          m_solver[idir]->m_hang       =  m_params.m_hang;
          m_solver[idir]->m_eps        =  m_params.m_tolerance;
          m_solver[idir]->m_verbosity  =  m_params.m_verbosity;
          m_solver[idir]->m_iterMax    =  m_params.m_iterMax;
          m_solver[idir]->m_normThresh =  m_params.m_normThresh;
          m_solver[idir]->setMGCycle(m_params.m_mgCycle);

          m_solver[idir]->init(phi,
                               rhs,
                               lmax,
                               lbase);

          if (m_params.m_orderTimeIntegration == 2)
            {
              m_tgaSolver[idir] = RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >
                (new AMRTGA<LevelData<EBCellFAB> >(m_solver[idir],
                                                   opFactory,
                                                   level0Dom,
                                                   m_params.m_refRatio,
                                                   numLevels,
                                                   m_params.m_verbosity));
            }
          else if (m_params.m_orderTimeIntegration == 1)
            {
              m_backwardSolver[idir] = RefCountedPtr<EBBackwardEuler>
                (new EBBackwardEuler(m_solver[idir],
                                     opFactory,
                                     level0Dom,
                                     m_params.m_refRatio,
                                     numLevels,
                                     m_params.m_verbosity));
            }
          else
            {
              MayDay::Error("EBAMRLinINS::defineProjections (2) -- bad order time integration");
            }

        }
      for (int ilev = 0; ilev < numLevels; ilev++)
        {
          delete phi[ilev];
          delete rhs[ilev];
        }
    }
}
/*********/
void EBAMRLinINS::
averageDown(Vector<LevelData<EBCellFAB>* >& a_data)
{
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }
      EBCoarseAverage& ebaverage = *avePtr;

      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } //end loop over levels

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      a_data[ilev]->exchange(interval);
    }
}
/*********/
void EBAMRLinINS::
averageDown(Vector<LevelData<EBFluxFAB>* >& a_data)
{
  //do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }
      EBCoarseAverage& ebaverage = *avePtr;

      ebaverage.average(*a_data[ilev-1],
                        *a_data[ilev  ],
                        interval);
    } //end loop over levels

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      a_data[ilev]->exchange(interval);
    }
}
/*********/
void EBAMRLinINS::
setCoveredStuffToZero(LevelData<EBCellFAB>& a_vort)
{
  for (DataIterator dit = a_vort.dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&  vortFAB =a_vort[dit()];
      Real covVal = 0.0;
      for (int icomp = 0; icomp < vortFAB.nComp(); icomp++)
        {
          vortFAB.setCoveredCellVal(covVal, icomp);
        }
    }
}
/*********/
void EBAMRLinINS::
defineIrregularData()
{
    for (int ilev = 0; ilev <=  m_finestLevel; ilev++)
    {
      m_coveredFaceLitLo[ilev]->define(m_grids[ilev]);
      m_coveredFaceLitHi[ilev]->define(m_grids[ilev]);
      m_coveredSetsLitLo[ilev]->define(m_grids[ilev]);
      m_coveredSetsLitHi[ilev]->define(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          Box litBox = m_grids[ilev].get(dit());
          litBox.grow(1);
          litBox &= m_domain[ilev];
          (*m_coveredFaceLitLo[ilev])[dit()].resize(SpaceDim);
          (*m_coveredFaceLitHi[ilev])[dit()].resize(SpaceDim);
          (*m_coveredSetsLitLo[ilev])[dit()].resize(SpaceDim);
          (*m_coveredSetsLitHi[ilev])[dit()].resize(SpaceDim);

          for (int idir = 0; idir < SpaceDim; idir++)
            {
              IntVectSet irregIVSPlus, irregIVSMinu;
              //get the covered sets and faces
              EBArith::computeCoveredFaces((*m_coveredFaceLitHi[ilev])[dit()][idir],
                                           (*m_coveredSetsLitHi[ilev])[dit()][idir],
                                           irregIVSPlus,idir, Side::Hi, m_ebisl[ilev][dit()], litBox);
              EBArith::computeCoveredFaces((*m_coveredFaceLitLo[ilev])[dit()][idir],
                                           (*m_coveredSetsLitLo[ilev])[dit()][idir],
                                           irregIVSMinu, idir, Side::Lo,  m_ebisl[ilev][dit()], litBox);

            }
        }

    }//end loop over levels

  for (int ilev = 0; ilev <=  m_finestLevel; ilev++)
    {
      m_coveredBaseAdvVelLo[ilev]->define(m_grids[ilev]);
      m_coveredBaseAdvVelHi[ilev]->define(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          (*m_coveredBaseAdvVelLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*m_coveredBaseAdvVelHi[ilev])[dit()].resize(SpaceDim, NULL);
          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              (*m_coveredBaseAdvVelLo[ilev])[dit()][idir]  = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, 1);
              (*m_coveredBaseAdvVelHi[ilev])[dit()][idir]  = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, 1);
            }
        }
    }
}
/*********/
void EBAMRLinINS::
postInitialize()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::postInitialize" << endl;
    }
}
/*********/
void EBAMRLinINS::
computeVorticity(LevelData<EBCellFAB>& a_vort, int a_level)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::computeVorticity" << endl;
    }
  CH_assert(m_isSetup);
  EBCellFactory ebcellfact(m_ebisl[a_level]);
  //define the data holder
  //vorticity is a vector in 3d, scalar in 2d
  int ncomp = 1;
  if (SpaceDim==3)
    {
      ncomp = SpaceDim;
    }
  a_vort.define(m_grids[a_level], ncomp, IntVect::Zero, ebcellfact);

  //interpolate velocity at coarse-fine interfaces
  //need some scratch space to do quadcfi

  int numLevels = m_finestLevel+1;
  Vector<LevelData<EBCellFAB>* > cellScratch;
  cellScratch.resize(numLevels, NULL);

// Sid: this is wrong if restarting from diff max level
/*
  for (int ilev=0; ilev < numLevels; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      cellScratch[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
    }
*/

  if (a_level > 0)
    {
// Sid added
      for (int ilev=0; ilev < numLevels; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          cellScratch[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
        }
// end add

      LevelData<EBCellFAB>& phiCoar = *cellScratch[a_level-1];
      LevelData<EBCellFAB>& phiFine = *cellScratch[a_level  ];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          Interval phiInterv(0, 0);
          Interval velInterv(idir, idir);
          m_velo[a_level-1]->copyTo(velInterv, phiCoar, phiInterv);
          m_velo[a_level  ]->copyTo(velInterv, phiFine, phiInterv);
          m_quadCFI[a_level]->interpolate(phiFine,
                                          phiCoar,
                                          phiInterv);

          //on copy back, we need ghost cells, so do the data iterator loop
          for (DataIterator dit = phiFine.dataIterator(); dit.ok(); ++dit)
            {
              Box region = phiFine[dit()].getRegion(); //includes ghost cells
              (*m_velo[a_level])[dit()].copy(region, velInterv, region, phiFine[dit()], phiInterv);
            }

          EBLevelDataOps::setVal(phiFine, 0.0);
          EBLevelDataOps::setVal(phiCoar, 0.0);
        }

    }
  m_velo[a_level]->exchange(Interval(0, SpaceDim-1));

  //clean up temp space
  for (int ilev=0; ilev < cellScratch.size(); ilev++)
    {
      if (cellScratch[ilev] != NULL)
        {
          delete cellScratch[ilev];
          cellScratch[ilev] = NULL;
        }
    }

  Real dxLev = m_dx[a_level];
  //do actual computation
  for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
    {
      EBCellFAB&  vortFAB =                a_vort[dit()];
      EBCellFAB&  veloFAB = (*m_velo[a_level])[dit()];
      const EBGraph& ebgraph =   m_ebisl[a_level][dit()].getEBGraph();
      vortFAB.setVal(0.);

      BaseFab<Real>& regVort = vortFAB.getSingleValuedFAB();
      BaseFab<Real>& regVelo = veloFAB.getSingleValuedFAB();
      Box interiorBox = m_grids[a_level].get(dit());
      interiorBox.grow(1);
      interiorBox &= m_domain[a_level];
      interiorBox.grow(-1);
      int vortIndex;
#if CH_SPACEDIM==3
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          vortIndex = idir;
#else
          vortIndex = 0;
          int idir = 3;
#endif

          //compute on all cells as regular
          FORT_COMPUTEVORT(CHF_FRA1(regVort, vortIndex),
                           CHF_CONST_FRA(regVelo),
                           CHF_BOX(interiorBox),
                           CHF_CONST_REAL(dxLev),
                           CHF_CONST_INT(idir));;

          //do cells on domain boundary as if they were irregular
          IntVectSet ivsIrreg(m_grids[a_level].get(dit()));
          ivsIrreg -= interiorBox;
          ivsIrreg |= ebgraph.getIrregCells(interiorBox);
          int diffDirVec[2];
          if (SpaceDim==2 || idir ==2)
            {
              diffDirVec[0] = 0;
              diffDirVec[1] = 1;
            }
          else if (idir == 0)
            {
              diffDirVec[0] = 1;
              diffDirVec[1] = 2;
            }
          else if (idir==1)
            {
              diffDirVec[0] = 2;
              diffDirVec[1] = 0;
            }
          else
            {
              MayDay::Error("missed a case");
            }

          for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
            {
              const VolIndex vof = vofit();
              Real vortValue = 0.0;

              //taking the derivative d(udvdir)/d(xdiffdir)
              for (int idiff = 0; idiff < 2; idiff++)
                {

                  Real signDiff = 0;
                  int diffDir= -1;
                  int velComp = -1;
                  if (idiff == 0)
                    {
                      signDiff = 1.0;
                      diffDir  = diffDirVec[0];
                      velComp  = diffDirVec[1];
                    }
                  else if (idiff == 1)
                    {
                      signDiff = -1.0;
                      diffDir  = diffDirVec[1];
                      velComp  = diffDirVec[0];
                    }
                  else
                    {
                      MayDay::Error("missed a case");
                    }

                  Vector<FaceIndex> hiFaces =  ebgraph.getFaces(vof, diffDir, Side::Hi);
                  Vector<FaceIndex> loFaces =  ebgraph.getFaces(vof, diffDir, Side::Lo);

                  bool hasHi = (hiFaces.size() == 1) && (!hiFaces[0].isBoundary());
                  bool hasLo = (loFaces.size() == 1) && (!loFaces[0].isBoundary());
                  Real diffValue = 0.0;
                  if (hasHi && hasLo)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =0.5*(hiValue - loValue)/dxLev;
                    }
                  else if (hasHi)
                    {
                      const VolIndex& vofHi = hiFaces[0].getVoF(Side::Hi);
                      Real hiValue = veloFAB(vofHi, velComp);
                      Real loValue = veloFAB(vof,   velComp);
                      diffValue =(hiValue - loValue)/dxLev;
                    }
                  else if (hasLo)
                    {
                      const VolIndex& vofLo = loFaces[0].getVoF(Side::Lo);
                      Real hiValue = veloFAB(vof  , velComp);
                      Real loValue = veloFAB(vofLo, velComp);
                      diffValue =  (hiValue - loValue)/dxLev;
                    }
                  else
                    {
                      diffValue = 0.0;
                    }
                  vortValue += signDiff*diffValue;
                } //end loop over idiff

              vortFAB(vof, vortIndex) = vortValue;

            } //end loop over irregular vofs
#if CH_SPACEDIM==3
        } //end loop over vort components in 3d
#endif

    }
}
/*********/
Real EBAMRLinINS::
computeDt()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::computeDt " << endl;
    }

  Real dt;
  if (m_useFixedDt)
    {
      dt = m_prescribedDt;
    }
  else
    {
      int numLevels = m_finestLevel+1;
      Real maxVel = 0.0;
      for (int ilev=0; ilev < numLevels; ilev++)
        {
          LevelData<EBCellFAB>& velNewLD = *m_velo[ilev];

          for (DataIterator dit = velNewLD.dataIterator(); dit.ok(); ++dit)
            {
              BaseFab<Real>& velFAB = (velNewLD[dit()]).getSingleValuedFAB();
              int iRegIrregCovered;
              CH_START(t1);
              const BaseFab<int>& maskFAB = m_ebisl[ilev][dit()].getEBGraph().getMask(iRegIrregCovered);
              CH_STOP(t1);
              if (iRegIrregCovered != -1)//not all covered
                {
                  if (iRegIrregCovered == 0)//has irreg
                    {
                      Box box = m_grids[ilev].get(dit());
                      for (int idir = 0; idir < SpaceDim; idir++)
                        {
                          FORT_MAXNORMMASK(CHF_REAL(maxVel),
                                           CHF_CONST_FRA1(velFAB,idir),
                                           CHF_BOX(box),
                                           CHF_CONST_FIA1(maskFAB,0));
                        }
                      IntVectSet ivs = m_ebisl[ilev][dit()].getMultiCells(box);
                      for (VoFIterator vofit(ivs, m_ebisl[ilev][dit()].getEBGraph()); vofit.ok(); ++vofit)
                        {
                          for (int idir = 0; idir < SpaceDim; idir++)
                            {
                              if (Abs(velNewLD[dit()](vofit(), idir)) > maxVel)
                              {
                                maxVel = Abs(velNewLD[dit()](vofit(), idir));
                              }
                            }
                        }
                    }
                  else//all reg
                    {
                      Box box = m_grids[ilev].get(dit());
                      for (int idir = 0; idir < SpaceDim; idir++)
                        {
                          FORT_MAXNORM(CHF_REAL(maxVel),
                                       CHF_CONST_FRA1(velFAB,idir),
                                       CHF_BOX(box));
                        }
                    }
                }
            }
        }
      CH_START(t2);
#ifdef CH_MPI
      Real tmp = 1.;
      int result = MPI_Allreduce(&maxVel, &tmp, 1, MPI_CH_REAL,
                                 MPI_MAX, Chombo_MPI::comm);
      if (result != MPI_SUCCESS)
        {
          MayDay::Error("communication error on norm");
        }
      maxVel = tmp;
#endif
      CH_STOP(t2);

      if (maxVel > 1.0e-10)
        {
          dt = m_params.m_cfl*m_dx[m_finestLevel]/maxVel;
        }
      else
        {
          dt = m_params.m_maxDt;
        }
    } //end if no prescribed dt

  //finally, enforce max dt grow factor
  if ( (m_dt > 0.) && (m_params.m_maxDtGrow > 0) &&
       (dt > m_params.m_maxDtGrow*m_dt) )
    {
      dt = m_params.m_maxDtGrow*m_dt;
    }

  return dt;
}
/*********/
Real EBAMRLinINS::
computeInitialDt()
{
  if (m_params.m_verbosity >= 3)
    {
      pout() << "EBAMRLinINS::computeInitialDt " << endl;
    }
  Real retval;
  if (m_useFixedDt)
    {
      retval = m_prescribedDt;
    }
  else
    {
      retval =  (m_params.m_initCFL/m_params.m_cfl)*computeDt();
    }
  return retval;
}
/*********/
void EBAMRLinINS::
extrapolateToCoveredFaces(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                          Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                          const Vector<LevelData<EBFluxFAB>* >&                 a_macOpen,
                          const Vector<LevelData<EBCellFAB>* >&                 a_cellOpen,
                          int                                                   a_idir,
                          Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd)
{
  if (a_ebLevAd == NULL)
  {
    a_ebLevAd = &m_ebLevAd;
  }
  // averageDown(a_macOpen);
  //extrapolate the projected velocity to covered faces
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      ProblemDomain curDomain(m_domain[ilev]);

      LevelData<EBFluxFAB>& faceValue = (LevelData<EBFluxFAB>&) *a_macOpen[ilev];
      faceValue.exchange(Interval(0,0));

      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {

          EBPatchAdvect& patcher  = (*a_ebLevAd)[ilev]->getPatchAdvect(dit());
          patcher.extrapToCoveredFaces(*(*a_coveredMacLo[ilev])[dit()][a_idir],
                                       (*a_macOpen[ilev])[dit()][a_idir],
                                       (*a_cellOpen[ilev])[dit()],
                                       (*m_coveredFaceLitLo[ilev])[dit()][a_idir],
                                       a_idir, Side::Lo, m_grids[ilev].get(dit()));

          patcher.extrapToCoveredFaces(*(*a_coveredMacHi[ilev])[dit()][a_idir],
                                       (*a_macOpen[ilev])[dit()][a_idir],
                                       (*a_cellOpen[ilev])[dit()],
                                       (*m_coveredFaceLitHi[ilev])[dit()][a_idir],
                                       a_idir, Side::Hi, m_grids[ilev].get(dit()));
        }
    }
}
/*********/
Real EBAMRLinINS::
run(Real a_maxTime, int a_maxStep)
{

}
/*********/
void EBAMRLinINS::
setupCoveredBaseAdvVelocity(const Vector<LevelData<EBCellFAB>* >& a_baseVelo,
                            const Vector<LevelData<EBFluxFAB>* >& a_baseAdvVelo)
{
  // just doing extrapolation to covered faces
  //    because the effect of source terms on extrap is hidden in a_adVelo
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    extrapolateToCoveredFaces(m_coveredBaseAdvVelLo,
                              m_coveredBaseAdvVelHi,
                              a_baseAdvVelo,
                              a_baseVelo, idir);
  }
}
/*********/
void EBAMRLinINS::
setupForStabilityRun(const Epetra_Vector&             a_x,
                     const Vector<DisjointBoxLayout>& a_baseflowDBL,
                     const Vector<EBLevelGrid>&       a_baseflowEBLG,
                     const std::string&               a_baseflowFile,
                     double                           a_pertScale,
                     bool                             a_incOverlapData,
                     bool                             a_setupForPlottingData)
{
  if (m_params.m_verbosity > 3)
  {
    pout() << "EBAMRLinINS::setupForStabilityRun" << endl;
  }

  // turn off regridding
  m_params.m_regridInterval = -1;
  m_params.m_doSFD = 0;
  m_isSetup = true;

  m_finestLevel = a_baseflowDBL.size() - 1;

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
    m_grids[ilev] = a_baseflowDBL[ilev];
  }

  defineEBISLs();
  defineNewVel();
  definePressure();
  defineProjections();

  // intialize data for StabilityRun:
  if (m_params.m_verbosity >= 3)
  {
    pout() << "EBAMRLinINS::initialize data for stability run" << endl;
  }

  EBAMRDataOps::setToZero(m_velo);
  EBAMRDataOps::setToZero(m_pres);
  EBAMRDataOps::setToZero(m_gphi);

  EBAMRDataOps::setToZero(m_baseVelo);
  EBAMRDataOps::setToZero(m_baseAdvVelo);

/*
  // Baseflow pressure grad
  Vector<LevelData<EBCellFAB>* > baseflowGPhi;
  baseflowGPhi.resize(m_finestLevel+1);

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
    baseflowGPhi[ilev] = new LevelData<EBCellFAB>();
    EBCellFactory ebcellfact(m_ebisl[ilev]);
    baseflowGPhi[ilev]->define(m_grids[ilev], SpaceDim,    IntVect::Zero, ebcellfact);
  }
*/
  // Do this only if !a_setupForPlottingData
#ifdef CH_USE_HDF5

  if (!a_setupForPlottingData)
  {
    HDF5Handle handleIn(a_baseflowFile, HDF5Handle::OPEN_RDONLY);
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      handleIn.setGroupToLevel(ilev);
      read<EBCellFAB>(handleIn, *m_baseVelo[ilev], "velo", m_grids[ilev], Interval(), false);
//      read<EBCellFAB>(handleIn, *baseflowGPhi[ilev], "gphi", m_grids[ilev], Interval(), false);
      read<EBFluxFAB>(handleIn, *m_baseAdvVelo[ilev], "advVel", m_grids[ilev], Interval(), false);
    }

    handleIn.close();
  }

#else

  MayDay::Error("EBAMRLinINS::setupForStabiltyRun needs HDF5 to read baseflowFile");

#endif

  // make U = a_pertScale*Uprime
  int nVeloComp = m_velo[0]->nComp();
//  int nPresComp = m_pres[0]->nComp();
//  int ntotComp = nVeloComp + nPresComp;
  int ntotComp = nVeloComp;

  ChomboEpetraOps::addEpetraVecToChomboData(m_velo, &a_x, 0., a_pertScale, 0, 0, nVeloComp, ntotComp, a_incOverlapData, m_params.m_refRatio);

//  ChomboEpetraOps::addEpetraVecToChomboData(m_pres, &a_x, 1., a_pertScale, 0, nVeloComp, nPresComp, ntotComp, a_incOverlapData, m_params.m_refRatio);

  // Average down finer levels onto coarser levels if a_incOverlapData is false
  // just averaging data onto coarse levels. The filtering and flux matching happens in run()
  if (!a_incOverlapData)
  {
    averageDown(m_velo);
    averageDown(m_pres);
  }

  defineIrregularData();

  if (!a_setupForPlottingData)
  {
    // setup baseflowAdvVelo if not setting up for plotting
    setupCoveredBaseAdvVelocity(m_baseVelo, m_baseAdvVelo);
  }

/*
  for (int ilev = 0; ilev <= baseflowGPhi.size(); ilev++)
  {
    delete baseflowGPhi[ilev];
    baseflowGPhi[ilev] = NULL;
  }
*/

  postInitialize();

//  m_doRestart = false; 
  m_doRestart = true;
  m_time = 0.;
}
/*********/
void EBAMRLinINS::
concludeStabilityRun(const std::string* a_pltName)
{

#ifdef CH_USE_HDF5

  writePlotFile(a_pltName);

#endif
}
/*********/
/*****************/
#ifdef CH_USE_HDF5
/*****************/
void
EBAMRLinINS::writePlotFile(const std::string* a_pltName)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::writePlotFile" << endl;
    }
  int curNumLevels = m_finestLevel + 1;

  Vector<string> presnames(1, string("pressure"));
#if CH_SPACEDIM == 2
  Vector<string> vortnames(1, string("vorticity"));
#else
  Vector<string> vortnames(SpaceDim);
#endif

  Vector<string> velonames(SpaceDim);
  Vector<string> gphinames(SpaceDim);
  Vector<string> names;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      char velochar[100];
      char gphichar[100];
      sprintf(velochar, "velocity%d", idir);
      sprintf(gphichar, "gradPres%d", idir);
      velonames[idir] = string(velochar);
      gphinames[idir] = string(gphichar);

#if CH_SPACEDIM==3
      char vortchar[100];
      sprintf(vortchar, "vorticity%d", idir);
      vortnames[idir] = string(vortchar);
#endif

    }

  names = velonames;
  names.append(gphinames);
  names.append(presnames);
  names.append(vortnames);

  int ncells = m_domain[0].size(0);

  bool replaceCovered = false;
  Vector<Real> coveredValues;

  int nlev = m_finestLevel + 1;
  Vector<LevelData<EBCellFAB>* > vorticity(nlev, NULL);
  Vector<LevelData<EBCellFAB>* > outputData(nlev, NULL);

  for (int ilev = 0; ilev < nlev; ilev++)
    {
      vorticity[ilev] = new LevelData<EBCellFAB>();
      computeVorticity(*vorticity[ilev], ilev);

#if CH_SPACEDIM==2
      int nvar = 2*SpaceDim + 2;
#else
      int nvar = 3*SpaceDim + 1;
#endif
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      outputData[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], nvar, IntVect::Zero, ebcellfact);

      Interval srcInterv, dstInterv;
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(0, SpaceDim-1);
      m_velo[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(SpaceDim, 2*SpaceDim-1);
      m_gphi[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim, 2*SpaceDim);
      m_pres[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

#if CH_SPACEDIM==2
      srcInterv = Interval(0, 0);
      dstInterv = Interval(2*SpaceDim+1, 2*SpaceDim+1);
      vorticity[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

#else
      srcInterv = Interval(0, SpaceDim-1);
      dstInterv = Interval(2*SpaceDim+1, 3*SpaceDim);
      vorticity[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);
#endif
      setCoveredStuffToZero(*outputData[ilev]);
    }

  string filename;
  if (a_pltName  == NULL)
  {
    filename = "plot.nx"+SSTR(ncells)+".step."+SSTR(m_curStep)+"."+SSTR(SpaceDim)+"d.hdf5";
  }
  else
  {
    filename = *a_pltName;
  }

  writeEBHDF5(filename,
              m_grids,
              outputData,
              names,
              m_domain[0].domainBox(),
              m_dx[0],
              m_dt,
              m_time,
              m_params.m_refRatio,
              curNumLevels,
              replaceCovered,
              coveredValues);

  for (int ilev = 0; ilev <nlev; ilev++)
    {
      delete vorticity[ilev];
      delete outputData[ilev];
    }
}
/*****************/
#endif //CH_USE_HDF5
/*****************/
