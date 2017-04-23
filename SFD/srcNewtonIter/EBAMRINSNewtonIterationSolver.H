/*
 *
 *
 *
 */

#ifndef _EBAMRINS_NEWTON_ITERATION_SOLVER_H_
#define _EBAMRINS_NEWTON_ITERATION_SOLVER_H_

#include "EBAMRNoSubcycle.H"
#include "EBAMRLinINS.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "DisjointBoxLayout.H"
#include "BiCGStabSolver.H"
#include "EBAMRINSNewtonIterationOp.H"

// class to perform Newton iterations

class EBAMRINSNewtonIterationSolver
{
public:

  EBAMRINSNewtonIterationSolver(const AMRParameters&               a_params,
                                const EBIBCFactory*                a_solverIBC,
                                const EBIBCFactory*                a_linIBC,
                                const ProblemDomain&               a_coarsestDomain,
                                Real                               a_nFlowSolverTime,
                                Real                               a_viscosity,
                                Real                               a_tol,
                                int                                a_maxIter,
                                const EBIndexSpace* const          a_ebisPtr = Chombo_EBIS::instance());

  EBAMRINSNewtonIterationSolver(const AMRParameters&               a_params,
                                const EBIBCFactory*                a_solverIBC,
                                const EBIBCFactory*                a_linIBC,
                                const Vector<DisjointBoxLayout>&   a_dbl,
                                const Vector<EBLevelGrid>&         a_eblg,
                                const Vector<Real>&                a_dx,
                                const ProblemDomain&               a_coarsestDomain,
                                Real                               a_nFlowSolverTime,
                                Real                               a_viscosity,
                                Real                               a_tol,
                                int                                a_maxIter,
                                const EBIndexSpace* const          a_ebisPtr = Chombo_EBIS::instance());

  virtual ~EBAMRINSNewtonIterationSolver();

  void oneStep(Vector<LevelData<EBCellFAB>* >&       a_stepChange,
               const Vector<LevelData<EBCellFAB>* >& a_baseState,
               const Vector<LevelData<EBFluxFAB>* >& a_baseAdvState,
               const Vector<LevelData<EBCellFAB>* >& a_baseStateChange);

  void solveNoInit(Vector<LevelData<EBCellFAB>* >&       a_finalState,
                   std::string                           a_restartFile,
                   Real                                  a_initialTime);

  void solveNoInit(Vector<LevelData<EBCellFAB>* >&       a_finalState,
                   const Vector<LevelData<EBCellFAB>* >& a_initialState);

  // number of initial steps to run
  //
  // 				WARNING!!!
  // 				WARNING!!! 
  // 				this calls a new on a_finalState
  // 				Need to call delete on a_finalState manually!!!!
  //
  // Try to call solveNoInit when possible !!!
  //
  void solve(Vector<LevelData<EBCellFAB>* >& a_finalState, int a_initSteps);

protected:

  AMRParameters          m_params;
  const EBIBCFactory*    m_ibcFact;
  const EBIBCFactory*    m_linIBCFact;

  Vector<DisjointBoxLayout> m_grids;
  Vector<EBLevelGrid>       m_eblg;
  const EBIndexSpace*       m_ebisPtr;
  Vector<Real>              m_dx;

  ProblemDomain             m_coarsestDomain;
  int  m_nlevels;
  Real m_viscosity;
  Real m_nFlowSolverTime;
  Real m_tol;
  int  m_maxIter;

  bool m_isSetupForNoInitRun;
  bool m_isSetup;

//  BiCGStabSolver<Vector<LevelData<EBCellFAB>* > > m_solver;

private:

  EBAMRINSNewtonIterationSolver()
  {;}
};

#endif
