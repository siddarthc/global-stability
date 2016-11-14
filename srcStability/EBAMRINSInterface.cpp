/*
 *
 *
 *
 */

#include "EBAMRINSInterface.H"

/*********/
/*********/
int EBAMRINSInterface::
getMaxLevelFromParmParse() const
{
  return m_params.m_maxLevel;
}
/*********/
/*********/
int EBAMRINSInterface::
getnEBGhost() const
{
  return 6; // to be consistent with EBAMRNoSubcycle::defineEBISLs()
}
/*********/
/*********/
int EBAMRINSInterface::
getnGhost() const
{
  return 3; // to be consistent with EBAMRNoSubcycle::defineNewVel()
}
/*********/
/*********/
int EBAMRINSInterface::
nComp() const
{
  int retval = SpaceDim + 1; // SpaceDim components of velocity + Pressure
}
/*********/
/*********/
#ifdef CH_USE_HDF5
/*********/
void EBAMRINSInterface::
readFileAndCopyToBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseflow, std::string a_baseflowFile, HDF5Handle& a_handleIn) const
{
  Interval veloInterv(0,SpaceDim-1);
  int nlevels = a_baseflow.size();

  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    a_handleIn.setGroupToLevel(ilev);
    ProblemDomain domain;
    getLevelDomain(&domain, ilev);
    const DisjointBoxLayout& dbl = a_baseflow[ilev]->getBoxes();
    EBLevelGrid eblg = EBLevelGrid(dbl, domain, getnEBGhost(), getEBISPtr());
    readFileAndCopyToBaseflow(a_baseflow[ilev], dbl, eblg.getEBISL(), domain, a_baseflowFile, a_handleIn);
  }

}
/*********/
/*********/
void EBAMRINSInterface::
readFileAndCopyToBaseflow(LevelData<EBCellFAB>* a_levBaseflow, const DisjointBoxLayout& a_levDBL, const EBISLayout& a_levEBISL, const ProblemDomain& a_domain, std::string a_baseflowFile, HDF5Handle& a_handleIn) const
{
  EBCellFactory ebcellfact(a_levEBISL);

  { // copy velocity
    int ncomp = SpaceDim;
    Interval srcInterv(0,SpaceDim-1);
    Interval dstInterv(0,SpaceDim-1);
    LevelData<EBCellFAB> tempLD(a_levDBL, ncomp, IntVect::Unit*getnGhost(), ebcellfact);
    read<EBCellFAB>(a_handleIn, tempLD, "velo", a_levDBL, Interval(), false);
    tempLD.copyTo(srcInterv, *a_levBaseflow, dstInterv);
  } // end velocity 

  { // copy pressure
    int ncomp = 1;
    Interval srcInterv(0,0);
    Interval dstInterv(SpaceDim, SpaceDim);
    LevelData<EBCellFAB> tempLD(a_levDBL, ncomp, IntVect::Unit*getnGhost(), ebcellfact);
    read<EBCellFAB>(a_handleIn, tempLD, "pres", a_levDBL, Interval(), false);
    tempLD.copyTo(srcInterv, *a_levBaseflow, dstInterv);
  }
}
/*********/
#endif
/*********/
/*********/
void EBAMRINSInterface::
computeSolution(Epetra_Vector& a_y, const Epetra_Vector& a_x, const Vector<LevelData<EBCellFAB>* >& a_baseflow, double a_eps, double a_integrationTime) const
{

}
