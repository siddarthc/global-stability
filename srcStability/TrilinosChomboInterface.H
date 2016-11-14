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
#include "Vector.H"
#include "EBCellFAB.H"
#include "LevelData.H"

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

  void setBaseflow(std::string a_baseflowFile);

  int nElementsOnThisProc() const;

  Epetra_Map getEpetraMap(const Epetra_Comm* a_commPtr) const;

  int computeL2Norm(const Epetra_Vector& a_v, double& a_result) const;

  int computeDotProd(const Epetra_Vector& a_v1, const Epetra_Vector& a_v2, double& a_result) const;

  void computeSolution(const Epetra_Vector& a_x, Epetra_Vector& a_y) const;

protected:

private:

  ChomboSolverInterface* m_solverInterface;

  // Baseflow parameters
  //

  Vector<LevelData<EBCellFAB>* > m_baseflow;


  // for weighted norms and inner products
  Epetra_Vector* m_volWeights;
  double         m_domainVolume;
  bool           m_isVolWeightsSet;

};

#endif
