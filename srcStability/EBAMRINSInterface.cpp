/*
 *
 *
 *
 */

#include "EBAMRINSInterface.H"
#include "ParmParse.H"

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

int EBAMRINSInterface::s_callCounter = -1;

/*********/
/*********/
EBAMRINSInterface::
EBAMRINSInterface(const AMRParameters& a_params,
                  const RefCountedPtr<EBIBCFactory>  a_IBC,
                  const ProblemDomain& a_coarsestDomain,
                  Real                 a_viscosity,
                  const EBIndexSpace* const a_ebisPtr)
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
/*********/
/*********/
int EBAMRINSInterface::
getMaxLevelFromParmParse() const
{
  CH_assert(m_isDefined);
  return m_params.m_maxLevel;
}
/*********/
/*********/
int EBAMRINSInterface::
getnEBGhost() const
{
  CH_assert(m_isDefined);
  return 6; // to be consistent with EBAMRNoSubcycle::defineEBISLs()
}
/*********/
/*********/
int EBAMRINSInterface::
getnGhost() const
{
  CH_assert(m_isDefined);
  return 3; // to be consistent with EBAMRNoSubcycle::defineNewVel()
}
/*********/
/*********/
int EBAMRINSInterface::
nComp() const
{
  CH_assert(m_isDefined);
  int retval = SpaceDim + 1; // SpaceDim components of velocity + Pressure
  return retval;
}
/*********/
/*********/
#ifdef CH_USE_HDF5
/*********/
void EBAMRINSInterface::
readFileAndCopyToBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseflow, std::string a_baseflowFile, HDF5Handle& a_handleIn) const
{
  CH_assert(m_isDefined);

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
  CH_assert(m_isDefined);

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
computeSolution(Epetra_Vector& a_y, const Epetra_Vector& a_x, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, const std::string& a_baseflowFile, double a_eps, double a_integrationTime, bool a_incOverlapData) const
{
  CH_assert(m_isDefined);

  int nlevels = a_baseflowDBL.size();
  int nComp = this->nComp();

  s_callCounter++;

  pout() << "doing y = EBAMRINSOp*x" << endl;

  // compute Frechet Derivative 
  
  // make a_yStar = (f(Ubar + eps*Uprime))/(2*eps)
  {
    EBAMRNoSubcycle solver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);
    solver.setupForStabilityRun(a_x, a_baseflowDBL, a_baseflowEBLG, a_baseflowFile, a_eps, a_incOverlapData);

    Real fixedDt = 0.;
    ParmParse pp;
    pp.query("fixed_dt", fixedDt);
    if (fixedDt > 1.e-12)
    {
      solver.useFixedDt(fixedDt);
    }

    int maxStep = 1000000;
    solver.run(a_integrationTime, maxStep);

    int check1 = a_y.PutScalar(0.);
    CH_assert(check1 == 0);

    const Vector<LevelData<EBCellFAB>* > veloSoln = solver.getVeloNew();
    const Vector<LevelData<EBCellFAB>* > presSoln = solver.getPresNew();

    int nVeloComp = veloSoln[0]->nComp();
    int nPresComp = presSoln[0]->nComp();
    int totComp = this->nComp();
    CH_assert(totComp == nVeloComp + nPresComp);

    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, veloSoln, 0., 1./(2.*a_eps), 0, 0, nVeloComp, totComp, a_incOverlapData, m_refRatio);
    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, presSoln, 0., 1./(2.*a_eps), nVeloComp, 0, nPresComp, totComp, a_incOverlapData, m_refRatio);

//    std::string pltName = "INS_soln_for_added_pert_at_step"+SSTR(s_callCounter)+".hdf5";
//    solver.concludeStabilityRun(&pltName);
  }

  // make a_y = a_yStar - (f(Ubar - eps*Uprime))/(2*eps)
  {
    EBAMRNoSubcycle solver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);
    solver.setupForStabilityRun(a_x, a_baseflowDBL, a_baseflowEBLG, a_baseflowFile, -1.*a_eps, a_incOverlapData);

    Real fixedDt = 0.;
    ParmParse pp;
    pp.query("fixed_dt", fixedDt);
    if (fixedDt > 1.e-12)
    {
      solver.useFixedDt(fixedDt);
    }

    int maxStep = 1000000;
    solver.run(a_integrationTime, maxStep);

    const Vector<LevelData<EBCellFAB>* > veloSoln = solver.getVeloNew();
    const Vector<LevelData<EBCellFAB>* > presSoln = solver.getPresNew();

    int nVeloComp = veloSoln[0]->nComp();
    int nPresComp = presSoln[0]->nComp();
    int totComp = this->nComp();
    CH_assert(totComp == nVeloComp + nPresComp);

    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, veloSoln, 1., -1./(2.*a_eps), 0, 0, nVeloComp, totComp, a_incOverlapData, m_refRatio);
    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, presSoln, 1., -1./(2.*a_eps), nVeloComp, 0, nPresComp, totComp, a_incOverlapData, m_refRatio);

//    std::string pltName = "INS_soln_for_subtracted_pert_at_step"+SSTR(s_callCounter)+".hdf5";
//    solver.concludeStabilityRun(&pltName);
  }

  pout() << "computed y = EMAMRINSOp*x" << endl;
  pout() << "Matrix-Vector Multiplication was computed " << s_callCounter << " times" << endl;
}
/*********/
/*********/
void EBAMRINSInterface::
plotEpetraVector(const Epetra_Vector& a_v, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_plotName, bool a_incOverlapData) const
{
  pout() << "plotting EigenEvector" << endl;
  EBAMRNoSubcycle solver(m_params, *m_ibcFact, m_coarsestDomain, m_viscosity);

  std::string baseflowFile = "dummyFile";

  solver.setupForStabilityRun(a_v, a_baseflowDBL, a_baseflowEBLG, baseflowFile, 1.0, a_incOverlapData, true);
  solver.run(0,0);
  solver.concludeStabilityRun(&a_plotName);
}
/*********/
/*********/
