/*
 *
 *
 *
 *
 */

#ifndef _EBAMRLININS_H_
#define _EBAMRLININS_H_

#include "EBAMRNoSubcycle.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBIBCFactory.H"
#include "AMRTGA.H"
#include "EBBackwardEuler.H"
#include "EBCompositeMACProjector.H"
#include "EBCompositeCCProjector.H"
#include "EBAMRPoissonOp.H"
#include "EBLevelAdvect.H"
#include "EBCoarseAverage.H"
#include "EBSimpleSolver.H"
#include "EBFastFR.H"
#include "ChomboSFDInterface.H"
#include "ChomboEpetraOps.H"

class EBAMRLinINS
{

public:

  EBAMRLinINS()
  {
    MayDay::Error("EBAMRLinINS::Invalid Constructor");
  }
  
  EBAMRLinINS(const AMRParameters&      a_params,
              const EBIBCFactory&       a_baseflowIBC,
              const EBIBCFactory&       a_solverIBC,
              const ProblemDomain&      a_coarsestDomain,
              Real                      a_viscosity,
              bool                      a_doAdjoint,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  virtual ~EBAMRLinINS();

  //fixes timestep to input
  void useFixedDt(Real a_dt);

  void setupForStabilityRun(const Epetra_Vector&             a_x, 
                            const Vector<DisjointBoxLayout>& a_baseflowDBL, 
                            const Vector<EBLevelGrid>&       a_baseflowEBLG, 
                            const std::string&               a_baseflowFile, 
                            double                           a_pertScale, 
                            bool                             a_incOverlapData, 
                            bool                             a_setupForPlottingData = false);

  /// advance solution.  returns time of final solution
  Real run(Real a_maxTime, int a_maxStep);

  Vector<LevelData<EBCellFAB>* >& getVeloNew()
  {
    return m_velo;
  }

  Vector<LevelData<EBCellFAB>* >& getPresNew()
  {
    return m_pres;
  }

  void concludeStabilityRun(const std::string* a_pltName);

protected:

  void allocateDataHolders();

  virtual void defineEBISLs();
  virtual void defineNewVel(const int a_startLevel = -1);
  virtual void definePressure(const int a_startLevel = -1);
  virtual void defineProjections();

  void averageDown(Vector<LevelData<EBCellFAB>* >&  a_data);
  void averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data);
  void setCoveredStuffToZero(LevelData<EBCellFAB>& a_vort);

  void postInitialize();

    /// compute dt
  virtual Real computeDt();

  /// compute dt with initial data
  Real computeInitialDt();

  // compute vorticity on level a_level
  void computeVorticity(LevelData<EBCellFAB>& a_vorticity, int a_level) ;

  void
  extrapolateToCoveredFaces(Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacLo,
                            Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >&  a_coveredMacHi,
                            const Vector<LevelData<EBFluxFAB>* >&                 a_macOpen,
                            const Vector<LevelData<EBCellFAB>* >&                 a_cellOpen,
                            int                                                   a_idir,
                            Vector<RefCountedPtr<EBLevelAdvect> >              *  a_ebLevAd = NULL);

#ifdef CH_USE_HDF5
  //
  virtual void writePlotFile(const std::string* a_pltName = NULL);

#endif

    //parmaters
  AMRParameters          m_params;
  EBIBC*                 m_baseflowIBC;
  EBIBC*                 m_solverIBC;
  const EBIndexSpace*    m_ebisPtr;

  //grid description
  Vector<ProblemDomain>     m_domain;
  Vector<Real>              m_dx;
  Vector<DisjointBoxLayout> m_grids;
  Vector<EBISLayout>        m_ebisl;
  Vector<EBLevelGrid>       m_eblg;

  EBCompositeCCProjector*        m_ccProjector;
  EBCompositeMACProjector*       m_macProjector;
  RefCountedPtr<AMRMultiGrid<LevelData<EBCellFAB> > >  m_solver[SpaceDim];
  BiCGStabSolver<LevelData<EBCellFAB> >          m_bottomSolver;
  EBSimpleSolver                                 m_bottomSolverSimp;
  RefCountedPtr<AMRTGA<LevelData<EBCellFAB> > >  m_tgaSolver[SpaceDim];
  RefCountedPtr<EBBackwardEuler>                 m_backwardSolver[SpaceDim];
  Vector<RefCountedPtr<EBQuadCFInterp> >    m_quadCFI;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveOper;
  Vector<RefCountedPtr<EBCoarseAverage> >   m_aveSpac;

  //  finest extant level
  int m_finestLevel;
  int m_curStep;

  //  time stepping information
  Real m_time;
  Real m_dt;
  Real m_prescribedDt;
  Real m_viscosity;

  long long m_pointsUpdated;

  //switches that i would rather lose but cannot
  bool m_useFixedDt;

  bool m_isSetup;
  //if viscosity == 0 this turns off viscous solvers, etc
  bool m_viscousCalc;

  bool m_doRestart;

  bool m_advanceGphiOnly;
  bool m_steadyState;
  bool m_stopAdvance;

  bool m_doAdjoint; 

  //data held in memory
  Vector<LevelData<EBCellFAB>* > m_velo;
  Vector<LevelData<EBCellFAB>* > m_pres;
  Vector<LevelData<EBCellFAB>* > m_gphi;

  Vector<RefCountedPtr<EBLevelAdvect>  > m_ebLevAd;
  Vector<RefCountedPtr<EBFastFR     >  > m_fluxReg;

  //baseflow data
  Vector<LevelData<EBCellFAB>* > m_baseVelo;
  Vector<LevelData<EBFluxFAB>* > m_baseAdvVelo; // baseflow advective velocity
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredBaseAdvVelLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredBaseAdvVelHi;
  
  //The following is the ugly part of this stuff.
  //Covered data (advection velocities, etc) are defined with this stuff.
  //These are used to define covered data and to iterate over it.
  //The actual data is defined on the fly to avoid confusion (no really)
  //and to save on memory since all these things are transitory.
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitLo;
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitHi;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitLo;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitHi;

  //defines the above stuff
  void defineIrregularData();

  void setupCoveredBaseAdvVelocity(const Vector<LevelData<EBCellFAB>* >& a_baseVelo,
                                   const Vector<LevelData<EBFluxFAB>* >& a_baseAdvVelo);

private:

  EBAMRLinINS(const EBAMRLinINS& a_input)
  {
    MayDay::Error("Invalid Constructor");
  }

  void operator=(const EBAMRLinINS& a_input)
  {
    MayDay::Error("Invalid Constructor");
  }

};

#endif 
