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
#include <Epetra_MpiComm.h>
//#include <Eperta_SerialComm.h>

class TrilinosSolverInterface
{
public:

  TrilinosSolverInterface()
  {
    m_isDtSizeSet = false;
    m_isEpsSet = false;
    m_isBaseflowSet = false;
  }

  virtual ~TrilinosSolverInterface()
  {

  }


  // save the steady state data (solution, grids, etc.)
  virtual void setBaseflow(std::string a_baseflowFile) = 0;

  virtual void setEps(double a_eps)
  {
    m_eps = a_eps;
    m_isEpsSet = true;
  }

  virtual void setIntegrationTime(double a_integrationTime)
  {
    m_integrationTime = a_integrationTime;
    m_isDtSizeSet = true;
  }

  // return the number of elements on this process
  // It helps in direct copying to and from epetra multivectors if the solver's map and the epetra map are one to one 
  virtual int nElementsOnThisProc() const = 0;

  // get the Epetra_Map
  //   the a_commPtr is not changed
  virtual Epetra_Map getEpetraMap(const Epetra_Comm* a_commPtr) const = 0;

  // compute L2 norm of the Epetra_MultiVec
  // return 0 is successfully computed
  virtual int computeL2Norm(const Epetra_Vector& a_v, double& a_result) const = 0;

  // compute the dot product of the MVs
  // return 0 is successfully computed
  virtual int computeDotProd(const Epetra_Vector& a_v1, const Epetra_Vector& a_v2, double& a_result) const = 0;

/*
  bool isLinearSolver()
  {
    return m_isLinearSolver;
  }
*/

  virtual void computeSolution(const Epetra_Vector& a_x, Epetra_Vector& a_y) const = 0;

  virtual bool isSetupForStabilityRun() const
  {
    return (m_isBaseflowSet && m_isEpsSet && m_isDtSizeSet);
  }

protected:
  
  double m_integrationTime;
  double m_eps;

  bool m_isLinearSolver;

  bool m_isBaseflowSet;
  bool m_isEpsSet;
  bool m_isDtSizeSet;
};

#endif
