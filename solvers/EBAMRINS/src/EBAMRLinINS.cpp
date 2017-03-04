/*
 *
 *
 *
 *
 */

#include "EBAMRLinINS.H"
#include "EBAMRDataOps.H"
#include "DirichletPoissonEBBC.H"

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

  m_baseVeloGrad.resize(SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    m_baseVeloGrad[idir].resize(nlevels, NULL);
  }

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

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      m_baseVeloGrad[idir][ilev] = new LevelData<EBCellFAB>();
    }
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

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      delete m_baseVeloGrad[idir][ilev];
    }
  }

  for (int ilev=0; ilev <= m_params.m_maxLevel; ilev++)
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
defineBaseflowVelocityGradients()
{
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
    EBCellFactory ebcellfact(m_ebisl[ilev]);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_baseVeloGrad[idir][ilev]->define(m_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact);
    }
  }

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    EBAMRDataOps::setToZero(m_baseVeloGrad[idir]);
  }
}
/*********/
void EBAMRLinINS::
computeBaseflowVelocityGradients()
{
  Vector<LevelData<EBCellFAB>* > velComp(m_finestLevel+1,NULL);
  
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
    EBCellFactory ebcellfact(m_ebisl[ilev]);
    velComp[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
  }

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      Interval srcInterv(idir, idir);
      Interval dstInterv(0, 0);
      m_baseVelo[ilev]->copyTo(srcInterv, *velComp[ilev], dstInterv);
    }

    m_ccProjector->gradient(m_baseVeloGrad[idir], velComp);
  } 

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
    delete velComp[ilev];
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
  defineBaseflowVelocityGradients();

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

  for (int idir = 0; idir < SpaceDim; idir++)
  {
    EBAMRDataOps::setToZero(m_baseVeloGrad[idir]);
  }

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
    computeBaseflowVelocityGradients();
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
void EBAMRLinINS::
computeExtraSourceForPredictor(Vector<LevelData<EBCellFAB>*> & a_source, const int& a_dir)
{
  
}
/*********/
void EBAMRLinINS::
computeExtraSourceForCorrector(Vector<LevelData<EBCellFAB>*> & a_source)
{

}
/*********/
void EBAMRLinINS::
transverseVelocityPredictor(Vector<LevelData<EBCellFAB>* >&    a_uDotDelU,
                            Vector<LevelData<EBCellFAB>* >&    a_scalOld,
                            bool                               a_reallyVelocity)
{
  int ncomp = a_uDotDelU[0]->nComp();
  int nlevels = m_finestLevel+1;

  //make temporaries with right number of variables
  Vector<LevelData<EBFluxFAB>* >                       macScratchVec(nlevels, NULL);
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScratchVecLo(nlevels, NULL);
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  coveredScratchVecHi(nlevels, NULL);

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      EBFluxFactory ebfluxfact(m_ebisl[ilev]);
      macScratchVec[ilev]       = new LevelData<EBFluxFAB>(m_grids[ilev], ncomp, 3*IntVect::Unit, ebfluxfact);

      coveredScratchVecLo[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      coveredScratchVecHi[ilev] = new LayoutData<Vector<BaseIVFAB<Real> * > >(m_grids[ilev]);
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {

          (*coveredScratchVecLo[ilev])[dit()].resize(SpaceDim, NULL);
          (*coveredScratchVecHi[ilev])[dit()].resize(SpaceDim, NULL);

          const EBGraph& ebgraph = m_ebisl[ilev][dit()].getEBGraph();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              (*coveredScratchVecLo[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitLo[ilev])[dit()][idir], ebgraph, ncomp);
              (*coveredScratchVecHi[ilev])[dit()][idir] = new BaseIVFAB<Real>((*m_coveredSetsLitHi[ilev])[dit()][idir], ebgraph, ncomp);
            }
        }
    }
  for (int icomp = 0; icomp < ncomp; icomp++)
    {
      EBPatchGodunov::setCurComp(icomp);
      EBPatchGodunov::setDoingAdvVel(0);

      RefCountedPtr<EBPhysIBCFactory> advectBC;
      if (a_reallyVelocity)
        {
          advectBC  = m_ibc->getVelAdvectBC(icomp);
          EBPatchGodunov::setDoingVel(1);
        }
      else
        {
          advectBC  = m_ibc->getScalarAdvectBC(icomp);
          EBPatchGodunov::setDoingVel(0);
        }

      for (int ilev=0; ilev <= m_finestLevel; ilev++)
        {
          Interval srcInterv(icomp, icomp);
          Interval dstInterv(0, 0);
          a_scalOld[ilev]->copyTo(srcInterv, *m_cellScratch[ilev], dstInterv);
        }

      //fill in viscous source term where necessary
      Vector<LevelData<EBCellFAB>* > * source = NULL;
      if (a_reallyVelocity && m_viscousCalc)
        {
          //cellscratch already holds velocity component

          DirichletPoissonEBBC::s_velComp = icomp;
          viscousSourceForAdvect(m_cellScratc2,   //returns holding source term = nu*lapl
                                 m_cellScratch,   //holds cell centered vel comp n
                                 m_cellScratc1,   //will hold the zero for the residual calc (zeroed inside routine)
                                 icomp,           //velocity component
                                 m_time);         //time, for BC setting
          source = &m_cellScratc2;
        }

      Vector<LevelData<EBCellFAB>*> extraSource;
      extraSource.resize(m_finestLevel+1, NULL);
      for (int ilev=0; ilev <= m_finestLevel; ilev++)
        {
          EBCellFactory ebcellfact(m_ebisl[ilev]);
          extraSource[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, 3*IntVect::Unit, ebcellfact);
        }
      EBAMRDataOps::setToZero(extraSource);

      computeExtraSourceForPredictor(extraSource, icomp);

      if (source != NULL)
        {
          for (int ilev=0; ilev<= m_finestLevel; ilev++)
            {
              for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
                {
                  (*(*source)[ilev])[dit()] += (*extraSource[ilev])[dit()];
                }
            }
        }
      else
        {
          source = &extraSource;
        }

      //cellscratch used for consstate
      extrapolateScalarCol(m_macScratch1,
                           m_coveredScratchLo,
                           m_coveredScratchHi,
                           advectBC,
                           m_advVel,
                           source,
                           m_cellScratch,
                           m_velo);

      for (int ilev=0; ilev<= m_finestLevel; ilev++)
        {
          delete extraSource[ilev];
        }

      if (a_reallyVelocity)
        {
          //correct with previously stored pressure gradient
          Real time = m_time + 0.5*m_dt;
          m_macProjector->setTime(time);
          m_macProjector->correctVelocityComponent(m_macScratch1, m_macGradient, icomp);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              int faceDir = idir;
              int velComp = icomp;
              //correct covered velocity with extrapolation of pressure gradient
              m_macProjector->correctVelocityComponent(m_coveredScratchLo,
                                                       m_coveredScratchHi,
                                                       m_coveredFaceLitLo,
                                                       m_coveredFaceLitHi,
                                                       m_coveredSetsLitLo,
                                                       m_coveredSetsLitHi,
                                                       m_macGradient, faceDir, velComp);
            }

          //overwrite extrapolated with advective velocity if appropriate
          //(the normal velocity here had the wrong boundary conditions)
          for (int ilev = 0; ilev <= m_finestLevel; ilev++)
            {
              for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
                {
                  Box box = grow(m_grids[ilev].get(dit()), 1);
                  Interval interv(0,0);
                  //copy over covered advective velocity
                  (*m_coveredScratchLo[ilev])[dit()][icomp]->copy(box, interv, box, *(*m_coveredAdvVelLo[ilev])[dit()][icomp], interv);
                  (*m_coveredScratchHi[ilev])[dit()][icomp]->copy(box, interv, box, *(*m_coveredAdvVelHi[ilev])[dit()][icomp], interv);

                  EBFluxFAB& macExtrapFAB    = (*m_macScratch1[ilev])[dit()];
                  const EBFluxFAB& advVelFAB = (*m_advVel[ilev])[dit()];
                  //icomp is the the component of the velocity and
                  //therefore also the face for which it is the normal component
                  macExtrapFAB[icomp].copy(advVelFAB[icomp]);
                }
            }
        }

      //copy extrapolated values into vector holders
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
        {
          Interval srcInterv(0, 0);
          Interval dstInterv(icomp, icomp);
          m_macScratch1[ilev]->copyTo(srcInterv, *macScratchVec[ilev],  dstInterv);
          int ibox = 0;
          for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
            {
              Box box = grow(m_grids[ilev].get(dit()), 1);
              box &= m_domain[ilev];
              for (int idir = 0; idir < SpaceDim; idir++)
                {
                  (*coveredScratchVecLo[ilev])[dit()][idir]->copy(box, dstInterv, box, *(*m_coveredScratchLo[ilev])[dit()][idir], srcInterv);
                  (*coveredScratchVecHi[ilev])[dit()][idir]->copy(box, dstInterv, box, *(*m_coveredScratchHi[ilev])[dit()][idir], srcInterv);
                }

              ibox++;
            }
        }

      //May need predicted velocity for EBAMRNoSubcycle extensions
      storePredictedVelocity(m_macScratch1,
                             m_coveredScratchLo,
                             m_coveredScratchHi,
                             icomp);

    } //end loop over velocity components  (icomp)

  //average down face centered stuff so that it makes sense at coarse-fine interfaces
  averageDown(macScratchVec);

  //compute the actual advective derivative
  //split this out to make it separately testable
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::computing advective derivative" << endl;
    }

  computeAdvectiveDerivative(a_uDotDelU, m_baseAdvVelo, macScratchVec,
                             m_coveredBaseAdvVelLo, m_coveredBaseAdvVelHi,
                             coveredScratchVecLo, coveredScratchVecHi);

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      delete macScratchVec[ilev];
      for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
        {
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              delete (*coveredScratchVecLo[ilev])[dit()][idir];
              delete (*coveredScratchVecHi[ilev])[dit()][idir];
            }
        }
      delete coveredScratchVecLo[ilev];
      delete coveredScratchVecHi[ilev];
    }
}
