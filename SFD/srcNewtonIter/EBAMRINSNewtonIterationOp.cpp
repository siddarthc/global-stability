/*
 *
 *
 *
 *
 */

#include "EBAMRINSNewtonIterationOp.H"
#include "EBAMRDataOps.H"
#include "GeneralUtils.H"

/*********/
EBAMRINSNewtonIterationOp::
EBAMRINSNewtonIterationOp(const AMRParameters&               a_params,
                          const Vector<LevelData<EBCellFAB>*>& a_baseVelo,
                          const Vector<LevelData<EBFluxFAB>*>& a_baseAdvVelo,
                          const EBIBCFactory*                a_solverIBC,
                          const Vector<DisjointBoxLayout>&   a_dbl,
                          const Vector<EBLevelGrid>&         a_eblg,
                          const Vector<Real>&                a_dx,
                          const ProblemDomain&               a_coarsestDomain,
                          Real                               a_nFlowSolverTime,
                          Real                               a_viscosity,
                          const EBIndexSpace* const          a_ebisPtr)
{
  m_coarsestDomain = a_coarsestDomain;
  m_ibcFact = a_solverIBC;
  m_grids = a_dbl;
  m_eblg = a_eblg;
  m_params = a_params;
  m_viscosity = a_viscosity;
  m_ebisPtr = a_ebisPtr;   
  m_nlevels = a_dbl.size();

  m_baseVelo = a_baseVelo;
  m_baseAdvVelo = a_baseAdvVelo;

  m_dx = a_dx;
  m_nFlowSolverTime = a_nFlowSolverTime;

  GeneralUtils::computeLinfNorm(m_basePertSize, m_baseVelo, m_dx, m_params.m_refRatio, false);
}
/**********/
EBAMRINSNewtonIterationOp::
~EBAMRINSNewtonIterationOp()
{

}
/*********/
void EBAMRINSNewtonIterationOp::
residual(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_phi, const Vector<LevelData<EBCellFAB>* >& a_rhs, bool a_homogeneous)
{

  EBAMRDataOps::setToZero(a_lhs);

  {
    EBAMRLinINS solver(m_params, *m_ibcFact, *m_ibcFact, m_coarsestDomain, m_viscosity, false);
    // normalize perturbation
    Real pertSize, scale;
    GeneralUtils::computeLinfNorm(pertSize, a_phi, m_dx, m_params.m_refRatio, false);
    scale = (pertSize > 1.e-3) ? m_basePertSize/pertSize : 1.;
    solver.setupForNewtonIterationRun(a_phi, m_baseVelo, m_baseAdvVelo, m_grids, m_eblg, scale);
    int maxSteps = 100000;
    solver.run(m_nFlowSolverTime, maxSteps);
    EBAMRDataOps::incr(a_lhs, solver.getVeloNew(), 1./scale);
    solver.concludeNewtonIterationRun();
  }

  EBAMRDataOps::incr(a_lhs, a_phi, -1.);
  EBAMRDataOps::incr(a_lhs, a_rhs, -1.);
}
/*********/
void EBAMRINSNewtonIterationOp::
preCond(Vector<LevelData<EBCellFAB>* >& a_cor, const Vector<LevelData<EBCellFAB>* >& a_residual)
{
  EBAMRDataOps::setToZero(a_cor);
  EBAMRDataOps::incr(a_cor, a_residual, 1.);
}
/*********/
void EBAMRINSNewtonIterationOp::
applyOp(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_phi, bool a_homogeneous)
{
  EBAMRDataOps::setToZero(a_lhs);

  {
    EBAMRLinINS solver(m_params, *m_ibcFact, *m_ibcFact, m_coarsestDomain, m_viscosity, false);
    // normalize perturbation
    Real pertSize, scale;
    GeneralUtils::computeLinfNorm(pertSize, a_phi, m_dx, m_params.m_refRatio, false);
    scale = (pertSize > 1.e-3) ? m_basePertSize/pertSize : 1.;
    solver.setupForNewtonIterationRun(a_phi, m_baseVelo, m_baseAdvVelo, m_grids, m_eblg, scale);
    int maxSteps = 100000;
    solver.run(m_nFlowSolverTime, maxSteps);
    EBAMRDataOps::incr(a_lhs, solver.getVeloNew(), 1./scale);
    solver.concludeNewtonIterationRun();
  }

  EBAMRDataOps::incr(a_lhs, a_phi, -1.); 
}
/*********/
void EBAMRINSNewtonIterationOp::
create(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_rhs)
{
//  CH_assert(a_lhs.size() == m_nlevels);

  a_lhs.resize(m_nlevels);
  int nComp = a_rhs[0]->nComp();
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    EBCellFactory ebcellfact(m_eblg[ilev].getEBISL());
    a_lhs[ilev] = new LevelData<EBCellFAB>();
    a_lhs[ilev]->define(m_eblg[ilev].getDBL(), nComp, a_rhs[ilev]->ghostVect(), ebcellfact);
  }
}
/*********/
void EBAMRINSNewtonIterationOp::
clear(Vector<LevelData<EBCellFAB>* >& a_lhs)
{
  for (int ilev = 0; ilev < m_nlevels; ilev++)
  {
    if (a_lhs[ilev] != NULL) delete a_lhs[ilev];
  }
}
/*********/
void EBAMRINSNewtonIterationOp::
assign(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_rhs)
{
  EBAMRDataOps::assign(a_lhs, a_rhs);
}
/*********/
Real EBAMRINSNewtonIterationOp::
dotProduct(const Vector<LevelData<EBCellFAB>* >& a_1, const Vector<LevelData<EBCellFAB>* >& a_2)
{
  Real retval;
  GeneralUtils::computeKappaDotProduct(retval, a_1, a_2, m_params.m_refRatio, m_dx[0]);
  return retval;
}
/*********/
void EBAMRINSNewtonIterationOp::
incr(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_x, Real a_scale)
{
  EBAMRDataOps::incr(a_lhs, a_x, a_scale);
}
/*********/
void EBAMRINSNewtonIterationOp::
axby(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_x, const Vector<LevelData<EBCellFAB>* >& a_y, Real a_a, Real a_b)
{
  EBAMRDataOps::axby(a_lhs, a_x, a_y, a_a, a_b);
}
/*********/
void EBAMRINSNewtonIterationOp::
scale(Vector<LevelData<EBCellFAB>* >& a_lhs, const Real& a_scale)
{
  EBAMRDataOps::scale(a_lhs, a_scale);
}
/*********/
Real EBAMRINSNewtonIterationOp::
norm(const Vector<LevelData<EBCellFAB>* >& a_rhs, int a_ord)
{
  Real retval;

  switch(a_ord)
  {
    case 0  : 
              GeneralUtils::computeLinfNorm(retval, a_rhs, m_dx, m_params.m_refRatio, false);
              break;
    case 1  :
              GeneralUtils::computeL1Norm(retval, a_rhs, m_dx, m_params.m_refRatio, false);
              break;
    case 2  :  
              GeneralUtils::computeL2Norm(retval, a_rhs, m_dx, m_params.m_refRatio, false);
              break;
    default :
              MayDay::Error("Invalid norm type");
  }

  return retval;
}
/*********/
void EBAMRINSNewtonIterationOp::
setToZero(Vector<LevelData<EBCellFAB>* >& a_lhs)
{
  EBAMRDataOps::setToZero(a_lhs);
}
/*********/
