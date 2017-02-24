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
              const EBIBCFactory&       a_IBC,
              const ProblemDomain&      a_coarsestDomain,
              Real                      a_viscosity,
              bool                      a_doAdjoint,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  virtual ~EBAMRLinINS();

protected:

  void allocateDataHolders();

    //parmaters
  AMRParameters          m_params;
  EBIBC*                 m_ibc;
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

  bool m_advanceGphiOnly;
  bool m_steadyState;
  bool m_stopAdvance;

  //data held in memory
  Vector<LevelData<EBCellFAB>* > m_velo;
  Vector<LevelData<EBCellFAB>* > m_pres;
  Vector<LevelData<EBCellFAB>* > m_gphi;

  Vector<RefCountedPtr<EBLevelAdvect>  > m_ebLevAd;
  Vector<RefCountedPtr<EBFastFR     >  > m_fluxReg;

  //The following is the ugly part of this stuff.
  //Covered data (advection velocities, etc) are defined with this stuff.
  //These are used to define covered data and to iterate over it.
  //The actual data is defined on the fly to avoid confusion (no really)
  //and to save on memory since all these things are transitory.
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitLo;
  Vector< LayoutData< Vector< Vector<VolIndex> > >* >  m_coveredFaceLitHi;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitLo;
  Vector< LayoutData< Vector< IntVectSet > >* >        m_coveredSetsLitHi;

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
