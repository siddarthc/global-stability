/*
 *
 *
 *
 */

#ifndef _EBAMRINS_NEWTON_ITERATION_OP_H_
#define _EBAMRINS_NEWTON_ITERATION_OP_H_

#include "EBAMRNoSubcycle.H"
#include "EBAMRLinINS.H"
#include "LinearSolver.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "DisjointBoxLayout.H"

// class to get the step change for Newton solver

class EBAMRINSNewtonIterationOp : public LinearOp< Vector< LevelData<EBCellFAB>* > >
{
public:

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
                            const EBIndexSpace* const          a_ebisPtr = Chombo_EBIS::instance());

  virtual ~EBAMRINSNewtonIterationOp();

  ///
  /**
     Say you are  solving L(phi) = rhs.   Make a_lhs = L(a_phi) - a_rhs.   If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void residual(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_phi, const Vector< LevelData<EBCellFAB>* >& a_rhs, bool a_homogeneous = false);

  ///
  /**
     Given the current state of the residual the correction, apply your preconditioner to a_cor.
   */
  virtual void preCond(Vector<LevelData<EBCellFAB>* >& a_cor, const Vector<LevelData<EBCellFAB>* >& a_residual);

  ///
  /**
     In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).  If a_homogeneous is true,
     evaluate the operator using homogeneous boundary conditions.
   */
  virtual void applyOp(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_phi, bool a_homogeneous = false);

  ///
  /**
     Creat data holder a_lhs that mirrors a_rhs.   You do not need to copy the data of a_rhs,
     just  make a holder the same size.
   */
  virtual void create(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_rhs);

  ///
  /**
     Opposite of create -- perform any operations required before lhs goes
     out of scope. In general, the only time this needs to be defined in
     a derived class is if the create() function called new. Otherwise, the
     default no-op function is sufficient.
  */
  virtual void clear(Vector<LevelData<EBCellFAB>* >& a_lhs);

  ///
  /**
     Set a_lhs  equal to a_rhs.
   */
  virtual void assign(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_rhs);

  ///
  /**
     Compute and return the dot product of a_1 and a_2.   In most contexts, this
     means return the sum over all data points of a_1*a_2.
   */
  virtual Real dotProduct(const Vector<LevelData<EBCellFAB>* >& a_1, const Vector<LevelData<EBCellFAB>* >& a_2);

  ///
  /**
     Increment by scaled amount (a_lhs += a_scale*a_x).
   */
  virtual void incr(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_x, Real a_scale);

  ///
  /**
     Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
   */
  virtual void axby(Vector<LevelData<EBCellFAB>* >& a_lhs, const Vector<LevelData<EBCellFAB>* >& a_x, const Vector<LevelData<EBCellFAB>* >& a_y, Real a_a, Real a_b);

  ///
  /**
     Multiply the input by a given scale (a_lhs *= a_scale).
   */
  virtual void scale(Vector<LevelData<EBCellFAB>* >& a_lhs, const Real& a_scale);

  ///
  /**
     Return the norm of  a_rhs.
     a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
   */
  virtual Real norm(const Vector<LevelData<EBCellFAB>* >& a_rhs, int a_ord);

  ///
  /**
     Set a_lhs to zero.
   */
  virtual void setToZero(Vector<LevelData<EBCellFAB>* >& a_lhs);

protected:

  AMRParameters          m_params;
  const EBIBCFactory*    m_ibcFact;
  
  Vector<LevelData<EBCellFAB>* > m_baseVelo;
  Vector<LevelData<EBFluxFAB>* > m_baseAdvVelo;

  Vector<DisjointBoxLayout> m_grids;
  Vector<EBLevelGrid>       m_eblg;
  const EBIndexSpace*       m_ebisPtr;
  Vector<Real>              m_dx;

  ProblemDomain             m_coarsestDomain;
  int  m_nlevels;
  Real m_viscosity;
  Real m_nFlowSolverTime;

  Real m_basePertSize;

private:

  EBAMRINSNewtonIterationOp()
  {;}
};

#endif
