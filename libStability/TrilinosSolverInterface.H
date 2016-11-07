/*
 *
 *
 *
 *
 */

#ifndef _TrilinosSolverInterface_H_
#define _TrilinosSolverInterface_H_

#include <string>
#include "Epetra_MultiVector.h"

class TrilinosSolverInterface
{
public:

  TrilinosSolverInterface()
  {

  }

  virtual ~TrilinosSolverInterface()
  {

  }


  // save the steady state data (solution, grids, etc.)
  virtual void setupForStabilityRun(std::string a_baseflowFile) = 0;

  virtual void setEps(double a_eps)
  {
    m_eps = a_eps;
  }

  virtual void setTimeStepSize(double a_timeStepSize)
  {
    m_timeStepSize = a_timeStepSize;
  }

  // return the number of elements on this process
  // It helps in direct copying to and from epetra multivectors is the solver's map and the epetra map are one to one 
  virtual int nElementsOnThisProc() const = 0;

  // compute L2 norm of the Epetra_MultiVec
  // return 0 is successfully computed
  virtual int computeL2Norm(const Epetra_MultiVector& a_mv, double* a_result) const = 0;

  // compute the dot product of the MVs
  // return 0 is successfully computed
  virtual int computeDotProd(const Epetra_MultiVector& a_mv1, const Epetra_MultiVector& a_mv2, double* a_result) const = 0;

  bool isLinearSolver()
  {
    return m_isLinearSolver;
  }

  virtual void computeSolution(const Epetra_MultiVector& a_x, Epetra_MultiVector& a_y) = 0;

protected:
  
  double m_timeStepSize;
  double m_eps;

  bool m_isLinearSolver;
};

#endif
