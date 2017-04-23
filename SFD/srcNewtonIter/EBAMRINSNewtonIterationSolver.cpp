/*
 *
 *
 *
 *
 */

#include "EBAMRINSNewtonIterationSolver.H"
#include "EBAMRDataOps.H"
#include "GeneralUtils.H"

/*********/
EBAMRINSNewtonIterationSolver::
EBAMRINSNewtonIterationSolver(const AMRParameters&               a_params,
                              const EBIBCFactory*                a_solverIBC,
                              const EBIBCFactory*                a_linIBC,
                              const ProblemDomain&               a_coarsestDomain,
                              Real                               a_nFlowSolverTime,
                              Real                               a_viscosity,
                              Real                               a_tol,
                              int                                a_maxIter,
                              const EBIndexSpace* const          a_ebisPtr)
{
  m_coarsestDomain = a_coarsestDomain;
  m_ibcFact = a_solverIBC;
  m_linIBCFact = a_linIBC;
  m_params = a_params;
  m_viscosity = a_viscosity;
  m_ebisPtr = a_ebisPtr;
  m_nFlowSolverTime = a_nFlowSolverTime;
  m_tol = a_tol;
  m_maxIter = a_maxIter;
  m_isSetup = true;
  m_isSetupForNoInitRun = false;
}
/*********/
EBAMRINSNewtonIterationSolver::
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
                              const EBIndexSpace* const          a_ebisPtr)
{
  m_coarsestDomain = a_coarsestDomain;
  m_ibcFact = a_solverIBC;
  m_linIBCFact = a_linIBC;
  m_grids = a_dbl;
  m_eblg = a_eblg;
  m_params = a_params;
  m_viscosity = a_viscosity;
  m_ebisPtr = a_ebisPtr;
  m_nlevels = a_dbl.size();
  m_dx = a_dx;
  m_nFlowSolverTime = a_nFlowSolverTime;
  m_tol = a_tol;
  m_maxIter = a_maxIter;
  m_isSetupForNoInitRun = true;
}
/*********/
EBAMRINSNewtonIterationSolver::
~EBAMRINSNewtonIterationSolver()
{

}
/*********/
void EBAMRINSNewtonIterationSolver::
oneStep(Vector<LevelData<EBCellFAB>* >&       a_stepChange,
        const Vector<LevelData<EBCellFAB>* >& a_baseState,
        const Vector<LevelData<EBFluxFAB>* >& a_baseAdvState,
        const Vector<LevelData<EBCellFAB>* >& a_baseStateChange)
{
  EBAMRINSNewtonIterationOp op(m_params, a_baseState, a_baseAdvState, m_linIBCFact, m_grids, m_eblg, m_dx, m_coarsestDomain, m_nFlowSolverTime, m_viscosity);

  GMRESSolver<Vector<LevelData<EBCellFAB>* > > solver;
  solver.define(&op, false);
  solver.m_verbosity = 5;
  solver.setRestartLen(64);
  solver.setConvergenceMetrics(0., m_tol);
  solver.solve(a_stepChange, a_baseStateChange);
}
/*********/
void EBAMRINSNewtonIterationSolver::
solveNoInit(Vector<LevelData<EBCellFAB>* >&       a_finalState,
            std::string                           a_restartFile,
            Real                                  a_initialTime)
{
  CH_assert(m_isSetupForNoInitRun);

  Real tol = 1.e6;
  int nIter = 0;

  Vector<LevelData<EBCellFAB>* > stepChange(m_nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > state(m_nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > stateChange(m_nlevels, NULL);
  Vector<LevelData<EBFluxFAB>* > advState(m_nlevels, NULL);

  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    stepChange[ilev] = new LevelData<EBCellFAB>();
    state[ilev] = new LevelData<EBCellFAB>();
    stateChange[ilev] = new LevelData<EBCellFAB>();
    advState[ilev] = new LevelData<EBFluxFAB>();

    EBCellFactory ebcellfact(m_eblg[ilev].getEBISL());
    stepChange[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
    state[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
    stateChange[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);

    EBFluxFactory ebfluxfact(m_eblg[ilev].getEBISL());
    advState[ilev]->define(m_grids[ilev], 1,  3*IntVect::Unit, ebfluxfact);
  }

  EBAMRDataOps::setToZero(stepChange);
  EBAMRDataOps::setToZero(state);
  EBAMRDataOps::setToZero(stateChange);
  EBAMRDataOps::setToZero(advState);

  {
    EBAMRNoSubcycle flowSolver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);
    flowSolver.setupForRestart(a_restartFile);
    EBAMRDataOps::incr(state, flowSolver.getVeloNew(), 1.);
    EBAMRDataOps::assign(advState, flowSolver.getAdvVelo());

    int maxSteps = 100000;
    flowSolver.run(a_initialTime + m_nFlowSolverTime, maxSteps);
    EBAMRDataOps::axby(stateChange, flowSolver.getVeloNew(), state, 1.0, -1.0);
    flowSolver.concludeNewtonIterationRun();
  }

  while ((tol > m_tol) &&  (nIter < m_maxIter))
  {  
    pout() << "doing Newton iteration no. " << nIter << endl;
    nIter++;

    EBAMRDataOps::setToZero(stepChange);

    oneStep(stepChange, state, advState, stateChange);
    EBAMRDataOps::incr(state, stepChange, -1.);
    GeneralUtils::computeL2Norm(tol, stepChange, m_dx, m_params.m_refRatio, false);

    if ((tol > m_tol) && (nIter < m_maxIter))
    {
      EBAMRNoSubcycle flowSolver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);
      flowSolver.setupForNewtonIterationRun(state, m_grids, m_eblg, false);

      // run for a small time to obtain the advection velocity
      int maxSteps = 100000;
      Real time = flowSolver.run(m_nFlowSolverTime, 2);
      EBAMRDataOps::setToZero(state);
      EBAMRDataOps::incr(state, flowSolver.getVeloNew(), 1.);
      EBAMRDataOps::assign(advState, flowSolver.getAdvVelo());

      flowSolver.run(m_nFlowSolverTime + time, maxSteps);
      EBAMRDataOps::axby(stateChange, flowSolver.getVeloNew(), state, 1.0, -1.0);
      flowSolver.concludeNewtonIterationRun();
    }  

    if (tol <= m_tol) 
    {
      pout() << "Newton iteration solver converged after " << nIter << " iterations" << endl;
    }
    else if (nIter >= m_maxIter)
    {
      pout() << "Newton iteration solver couldn't find root after " << nIter << " iterations" << endl;
      pout() << "Newton iteration solver L2 norm of step change is " << tol << endl;
    }
  }

  EBAMRDataOps::setToZero(a_finalState);
  EBAMRDataOps::incr(a_finalState, state, 1.);

  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    if (stepChange[ilev] != NULL) delete stepChange[ilev];
    if (state[ilev] != NULL) delete state[ilev];
    if (stateChange[ilev] != NULL) delete stateChange[ilev];
    if (advState[ilev] != NULL) delete advState[ilev];
  }
}
/*********/
void EBAMRINSNewtonIterationSolver::
solveNoInit(Vector<LevelData<EBCellFAB>* >&       a_finalState,
            const Vector<LevelData<EBCellFAB>* >& a_initialState)
{
  CH_assert(m_isSetupForNoInitRun);

  Real tol = 1.e6;
  int nIter = 0;

  Vector<LevelData<EBCellFAB>* > stepChange(m_nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > state(m_nlevels, NULL);
  Vector<LevelData<EBCellFAB>* > stateChange(m_nlevels, NULL);
  Vector<LevelData<EBFluxFAB>* > advState(m_nlevels, NULL);

  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    stepChange[ilev] = new LevelData<EBCellFAB>();
    state[ilev] = new LevelData<EBCellFAB>();
    stateChange[ilev] = new LevelData<EBCellFAB>();
    advState[ilev] = new LevelData<EBFluxFAB>();

    EBCellFactory ebcellfact(m_eblg[ilev].getEBISL());
    stepChange[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
    state[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
    stateChange[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);

    EBFluxFactory ebfluxfact(m_eblg[ilev].getEBISL());
    advState[ilev]->define(m_grids[ilev], 1,  3*IntVect::Unit, ebfluxfact);
  }

  EBAMRDataOps::setToZero(stepChange);
  EBAMRDataOps::setToZero(state);
  EBAMRDataOps::setToZero(stateChange);
  EBAMRDataOps::setToZero(advState);

  EBAMRDataOps::incr(state, a_initialState, 1.);

  while ((tol > m_tol) &&  (nIter < m_maxIter))
  {  
    pout() << "doing Newton iteration no. " << nIter << endl;
    nIter++;

    EBAMRDataOps::setToZero(stateChange);
    EBAMRDataOps::setToZero(stepChange);

    {
      EBAMRNoSubcycle flowSolver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);
      flowSolver.setupForNewtonIterationRun(state, m_grids, m_eblg, false);
    
      // run for a small time to obtain the advection velocity
      int maxSteps = 100000;
      Real time = flowSolver.run(m_nFlowSolverTime, 2);
      EBAMRDataOps::setToZero(state);
      EBAMRDataOps::incr(state, flowSolver.getVeloNew(), 1.);
      EBAMRDataOps::assign(advState, flowSolver.getAdvVelo());

      flowSolver.run(m_nFlowSolverTime, maxSteps);
      EBAMRDataOps::axby(stateChange, flowSolver.getVeloNew(), state, 1.0, -1.0);
      flowSolver.concludeNewtonIterationRun();
    }

    oneStep(stepChange, state, advState, stateChange);
    EBAMRDataOps::incr(state, stepChange, -1.);
    GeneralUtils::computeL2Norm(tol, stepChange, m_dx, m_params.m_refRatio, false);

    if (tol < m_tol) 
    {
      pout() << "Newton iteration solver converged after " << nIter << " iterations" << endl;
    }
    else if (nIter >= m_maxIter)
    {
      pout() << "Newton iteration solver couldn't find root after " << nIter << " iterations" << endl;
      pout() << "Newton iteration solver L2 norm of step change is " << tol << endl;
    }
  }

  EBAMRDataOps::setToZero(a_finalState);
  EBAMRDataOps::incr(a_finalState, state, 1.);

  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    if (stepChange[ilev] != NULL) delete stepChange[ilev];
    if (state[ilev] != NULL) delete state[ilev];
    if (stateChange[ilev] != NULL) delete stateChange[ilev];
    if (advState[ilev] != NULL) delete advState[ilev];
  }
}
/*********/
void EBAMRINSNewtonIterationSolver::
solve(Vector<LevelData<EBCellFAB>* >& a_finalState, int a_initSteps)
{
  CH_assert(m_isSetup);

  m_nlevels = m_params.m_maxLevel+1;
  m_grids.resize(m_nlevels);
  m_eblg.resize(m_nlevels);
  m_dx.resize(m_nlevels);
  Vector<LevelData<EBCellFAB>* > initialState(m_nlevels, NULL);

  a_finalState.resize(m_nlevels);

  // run to get initial stuff
  {
    EBAMRNoSubcycle flowSolver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);
    flowSolver.setupForAMRRun();
    Real maxTime = 100000.0;
    flowSolver.run(maxTime, a_initSteps);

    m_grids = flowSolver.getGrids();
    m_eblg = flowSolver.getEBLG();
    m_dx = flowSolver.getDx();

    for (int ilev = 0; ilev < m_nlevels; ilev++)
    {
      initialState[ilev] = new LevelData<EBCellFAB>();
      a_finalState[ilev] = new LevelData<EBCellFAB>();  // not safe!!!
      EBCellFactory ebcellfact(m_eblg[ilev].getEBISL());
      initialState[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
      a_finalState[ilev]->define(m_grids[ilev], SpaceDim, 3*IntVect::Unit, ebcellfact);
    }

    EBAMRDataOps::setToZero(initialState);
    EBAMRDataOps::setToZero(a_finalState);
    EBAMRDataOps::incr(initialState, flowSolver.getVeloNew(), 1.0);
  }

  m_isSetupForNoInitRun = true;
  solveNoInit(a_finalState, initialState);

  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    if (initialState[ilev] != NULL) delete initialState[ilev];
  }
}
