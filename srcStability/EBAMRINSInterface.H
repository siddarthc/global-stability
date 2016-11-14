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

  // need AMRINSParams, EBISPtr, coarsestDomain, EBISPtr, etc.
  EBAMRINSInterface()
  {

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

#endif

  virtual void computeSolution(Epetra_Vector& a_y, const Epetra_Vector& a_x, const Vector<LevelData<EBCellFAB>* >& a_baseflow, double a_eps, double a_integrationTime) const;

  virtual void readFileAndCopyToBaseflow(LevelData<EBCellFAB>* a_levBaseflow, const DisjointBoxLayout& a_levDBL, const EBISLayout& a_levEBISL, const ProblemDomain& a_domain, std::string a_baseflowFile, HDF5Handle& a_handleIn) const;

private:

  AMRParameters m_params;
};
#endif
