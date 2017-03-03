/*
 *
 *
 *
 *
 */

#include "EBAMRLinINS.H"
#include "EBAMRDataOps.H"

/*********/
EBAMRLinINS::
EBAMRLinINS(const AMRParameters&      a_params,
            const EBIBCFactory&       a_baseflowIBCFact,
            const EBIBCFactory&       a_solverIBCFact,
            const ProblemDomain&      a_coarsestDomain,
            Real                      a_viscosity,
            bool                      a_doAdjoint,
            const EBIndexSpace* const a_ebisPtr) : 
EBAMRNoSubcycle(a_params, a_solverIBCFact, a_coarsestDomain, a_viscosity, a_ebisPtr)
{
    if (a_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::EBAMRLinINS" << endl;
    }

  m_doAdjoint = a_doAdjoint;

  //create initial and boundary condition object
  m_baseflowIBC    = a_baseflowIBCFact.create();

  int nlevels = m_params.m_maxLevel + 1;

    // base flow data
  m_baseVelo.resize(nlevels, NULL);
  m_baseAdvVelo.resize(nlevels, NULL);
  m_coveredBaseAdvVelLo.resize(nlevels, NULL);
  m_coveredBaseAdvVelHi.resize(nlevels, NULL);

  allocateBaseflowDataHolders();  
}
/*********/
void EBAMRLinINS::
allocateBaseflowDataHolders()
{
  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
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

  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
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
}
/*********/
void EBAMRLinINS::
defineNewBaseflowVel(const int a_startLevel)
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::defineNewBaseflowVel" << endl;
    }

  int startLevel = 0;
  if (a_startLevel > startLevel)
    {
      startLevel=a_startLevel;
    }
  for (int ilev = startLevel; ilev <= m_finestLevel; ilev++)
    {
      EBCellFactory ebcellfact(m_ebisl[ilev]);
      m_baseVelo[ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      m_baseAdvVelo[ilev]->define(m_grids[ilev], 1,  3*IntVect::Unit, ebfluxfact);
    }
  EBAMRDataOps::setToZero(m_baseVelo);
  EBAMRDataOps::setToZero(m_baseAdvVelo);
}
/*********/
void EBAMRLinINS::
defineBaseflowIrregularData()
{
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

  defineNewBaseflowVel();

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

  defineBaseflowIrregularData();

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
conclude()
{

}
/*********/
