/*
 *
 *
 *
 */

#ifndef _Chombo_Solver_Interface_H_
#define _Chombo_Solver_Interface_H_

#include "EBIndexSpace.H"

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
  virtual void setLevelDomain(ProblemDomain* a_domain, int a_level) const
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

  // read from the file a_baseflowFile and copy to a_baseflow
  // 	a_baseflow need not be defined in this routine
  virtual void readFileAndCopyToBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseflow, std::string& a_baseflowFile) const = 0;
  
protected:

  ProblemDomain       m_coarsestDomain;
  Real                m_coarsestDx;
  Vector<Real>        m_refRatio;
  const EBIndexSpace* m_ebisPtr; 

private:

};

#endif
