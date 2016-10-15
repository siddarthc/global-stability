/*
 *
 *
 *
 *
 */

#ifndef _StabilityEvaluator_H_
#define _StabilityEvaluator_H_

#include "SolverFactory.H"

using namespace std;

template <typename T>
class StabilityEvaluator
{

public:

  StabilityEvaluator()
  {
    m_isDefined = false;
  }

  StabilityEvaluator(SolverFactory& a_solverFact, int a_subspaceDim, double a_eps, double a_timeStep)
  {
    define(a_solver, a_baseFlow, a_subspaceDim, a_eps, a_timeStep);
  } 

  void define(SolverFactory& a_solverFact, int a_subspaceDim, double a_eps, double a_timeStep)
  {
    m_solver = a_solverFact.create();
    
    m_eps         = a_eps;
    m_subspaceDim = a_subspaceDim;
    m_timeStep    = a_timeStep;

    m_isDefined = true;
  }

  ~StabilityEvaluator()
  {

  }

private:
  
  bool m_isDefined;

  Solver* m_solver;

  int     m_subspaceDim;
  double  m_eps;
  double  m_timeStep;
};

#endif
