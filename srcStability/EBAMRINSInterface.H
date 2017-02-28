/*
 *
 *
 *
 */

#ifndef _EBAMRINS_Interface_H_
#define _EBAMRINS_Interface_H_

#include "ChomboSolverInterface.H"
#include "EBAMRNoSubcycle.H"

class EBAMRINSInterface : public ChomboSolverInterface
{
public:

  // a_incLevOverlapData is a flag to determine how data is copied to and from Epetra objects
  //   if a_incLevelOverlapData = true, it includes the overlap data
  //   if false, it only writes fines level's data in the overlap region
  // a_doFirstOrderDeriv is a flag for order of Frechet Derivative
  EBAMRINSInterface(const AMRParameters& a_params,
                    const RefCountedPtr<EBIBCFactory>& a_baseflowIBC,
                    const RefCountedPtr<EBIBCFactory>& a_solverIBC, 
                    const ProblemDomain& a_coarsestDomain,
                    Real                 a_viscosity,
                    bool                 a_plotSnapshots,
                    bool                 a_doLinINS,
                    bool                 a_doFirstOrderFreDeriv = false,
                    const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());

  virtual ~EBAMRINSInterface()
  {

  }

  virtual int getMaxLevelFromParmParse() const;

  virtual int getnEBGhost() const;

  virtual int getnGhost() const;

  virtual int nComp() const;

#ifdef CH_USE_HDF5

  virtual void readFileAndCopyToBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseflow, std::string a_baseflowFile, HDF5Handle& a_handleIn) const;

  virtual void readFileAndCopyToBaseflow(LevelData<EBCellFAB>* a_levBaseflow, const DisjointBoxLayout& a_levDBL, const EBISLayout& a_levEBISL, const ProblemDomain& a_domain, std::string a_baseflowFile, HDF5Handle& a_handleIn) const;
#endif

  virtual void computeSolution(Epetra_Vector& a_y, const Epetra_Vector& a_x, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, const std::string& a_baseflowFile, double a_eps, double a_integrationTime, bool a_incOverlapData) const;

  virtual void plotEpetraVector(const Epetra_Vector& a_v, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_plotName, bool a_incOverlapData) const;

  virtual void getBaseflow(Epetra_Vector& a_v, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_baseflowFile, bool a_incOverlapData) const;

  void getBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseLD, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_baseflowFile, bool a_incOverlapData) const;

private:

  AMRParameters m_params;
  RefCountedPtr<EBIBCFactory> m_baseflowIBCFact;
  RefCountedPtr<EBIBCFactory> m_solverIBCFact;
  Real m_viscosity;
  bool m_isDefined;
  bool m_plotSnapshots;
  bool m_doLinINS;
  bool m_doFirstOrderFreDeriv;

  Real m_baseflowL2Norm;

  static int s_callCounter; // a counter for no. of times computeSolution is called

  EBAMRINSInterface()
  {
    MayDay::Error("invalid constructor");
  }

  EBAMRINSInterface(const EBAMRINSInterface& a_input)
  {
    MayDay::Error("invalid constructor");
  }

  void operator=(const EBAMRINSInterface& a_input)
  {
    MayDay::Error("invalid constructor");
  }
};
#endif
