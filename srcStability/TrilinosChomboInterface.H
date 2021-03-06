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
#include "EBLevelGrid.H"
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
  
  // if a_incOverlapData is true, the data in coarse-fine overlap region into Epetra Vector is also copied
  // if a_doVolWeighting is true, volume-weighted norms and inner-products are used
  TrilinosChomboInterface(const RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact, bool a_incOverlapData, bool a_doWeighting);

  void setBaseflow(const std::string& a_baseflowFile, const Epetra_Comm* a_commPtr);

  int nElementsOnThisProc() const;

  Epetra_Map getEpetraMap(const Epetra_Comm* a_commPtr) const;

  int computeL2Norm(const Epetra_Vector& a_v, double& a_result) const;

  int computeDotProd(const Epetra_Vector& a_v1, const Epetra_Vector& a_v2, double& a_result) const;

  void computeSolution(const Epetra_Vector& a_x, Epetra_Vector& a_y) const;

  void plotEpetraVector(const Epetra_Vector& a_v, std::string a_plotName) const;

protected:

private:

  void computeBaseflowNorm(const Epetra_Comm* a_commPtr);

  ChomboSolverInterface* m_solverInterface;

  bool                   m_incOverlapData;
  bool                   m_doWeighting;
  // Baseflow parameters
  //
  std::string               m_baseflowFile;
  Vector<DisjointBoxLayout> m_baseflowDBL;
  Vector<EBLevelGrid> m_baseflowEBLG;

  // for weighted norms and inner products
  Epetra_Vector* m_volWeights;
  double         m_domainVolume;
  bool           m_isVolWeightsSet;

  // Baseflow L2 Norm
  bool m_isBaseflowNormSet;
  Real m_baseflowL2Norm;
};

#endif
