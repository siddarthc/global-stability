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
  EBAMRINSInterface(const AMRParameters& a_params,
                    const RefCountedPtr<EBIBCFactory>  a_IBC,
                    const ProblemDomain& a_coarsestDomain,
                    Real                 a_viscosity,
                    const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance())
  {
    m_isDefined = true;
    m_params = a_params;
    m_ibcFact = a_IBC;
    m_coarsestDomain = a_coarsestDomain;
    m_viscosity = a_viscosity;
    m_ebisPtr = a_ebisPtr;
    m_refRatio = m_params.m_refRatio;
    m_coarsestDx = m_params.m_domainLength/Real(a_coarsestDomain.size(0));
  }

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

private:

  AMRParameters m_params;
  RefCountedPtr<EBIBCFactory>  m_ibcFact;
  Real m_viscosity;
  bool m_isDefined;

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
