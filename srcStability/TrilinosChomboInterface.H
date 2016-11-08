/*
 *
 *
 *
 *
 */

#ifndef _TRILINOS_CHOMBO_INTERFACE_H_
#define _TRILINOS_CHOMBO_INTERFACE_H_

#include "TrilinosSolverInterface.H"
#include "MayDay.H"
#include "ChomboSolverInterfaceFactory.H"
#include "RefCountedPtr.H"

class TrilinosChomboInterface : public TrilinosSolverInterface
{
public:

  // bad construction
  TrilinosChomboInterface()
  {
    MayDay::Error("TrilinosChomboInterface::invalid constructor");
  }

  virtual ~TrilinosChomboInterface();
  
  TrilinosChomboInterface(const RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact);

  void setupForStabilityRun(std::string a_baseflowFile);

  int nElementsOnThisProc() const;

  int computeL2Norm(const Epetra_MultiVector& a_mv, double* a_result) const;

  int computeDotProd(const Epetra_MultiVector& a_mv1, const Epetra_MultiVector& a_mv2, double* a_result) const;

  void computeSolution(const Epetra_MultiVector& a_x, Epetra_MultiVector& a_y) const;

protected:

private:

  ChomboSolverInterface* m_solverInterface;

};

#endif
