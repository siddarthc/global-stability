/*
 *
 *
 *
 */

#ifndef _Chombo_Solver_Interface_H_
#define _Chombo_Solver_Interface_H_

#include "EBIndexSpace.H"
#include <Epetra_Vector.h>

class ChomboSolverInterface
{
public:

  ChomboSolverInterface()
  {  

  }

  virtual ~ChomboSolverInterface()
  {

  }

  // return the maximum level from the input file
  virtual int getMaxLevelFromParmParse() const = 0;

  // set the ProblemDomain for the given level
  virtual void getLevelDomain(ProblemDomain* a_domain, int a_level) const
  {
    if (a_level > 0)
    {
      CH_assert(m_refRatio.size() >= a_level-1);
    }

    *a_domain = m_coarsestDomain;
    for (int ilev = 1; ilev <= a_level; ilev++)
    {
      a_domain->refine(m_refRatio[ilev-1]);
    } 
  }

  // get EBIS*
  virtual const EBIndexSpace* getEBISPtr() const
  {
    CH_assert(m_ebisPtr != NULL);
    return m_ebisPtr;
  }

  // get numEBGhost: see defineEBISLs() in EBAMRNoSubcycle
  virtual int getnEBGhost() const = 0;

  // get the number of ghost cells in defining the baseflow data
  virtual int getnGhost() const = 0;

  // get nComp in baseflow vector
  virtual int nComp() const = 0;

#ifdef CH_USE_HDF5

  // read from the file a_baseflowFile and copy to a_baseflow
  // 	a_baseflow need not be defined in this routine
  virtual void readFileAndCopyToBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseflow, std::string a_baseflowFile, HDF5Handle& a_handleIn) const = 0;

  virtual void readFileAndCopyToBaseflow(LevelData<EBCellFAB>* a_levBaseflow, const DisjointBoxLayout& a_levDBL, const EBISLayout& a_levEBISL, const ProblemDomain& a_domain, std::string a_baseflowFile, HDF5Handle& a_handleIn) const = 0;

#endif

  // compute the solution (do: a_y = A*a_x)
  virtual void computeSolution(Epetra_Vector& a_y, const Epetra_Vector& a_x, const Vector<LevelData<EBCellFAB>* >& a_baseflow, double a_eps, double a_integrationTime) const = 0;
  
protected:

  ProblemDomain       m_coarsestDomain;
  Real                m_coarsestDx;
  Vector<Real>        m_refRatio;
  const EBIndexSpace* m_ebisPtr; 

private:

};

#endif
