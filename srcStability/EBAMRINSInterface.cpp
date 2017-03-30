/*
 *
 *
 *
 */

#include "EBAMRINSInterface.H"
#include "ParmParse.H"
#include "EBAMRDataOps.H"

#include "EBAMRLinINS.H"

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

int EBAMRINSInterface::s_callCounter = -1;

/*********/
/*********/
EBAMRINSInterface::
EBAMRINSInterface(const AMRParameters& a_params,
                  const RefCountedPtr<EBIBCFactory>& a_baseflowIBC,
                  const RefCountedPtr<EBIBCFactory>& a_solverIBC,
                  const ProblemDomain& a_coarsestDomain,
                  Real                 a_viscosity,
                  bool                 a_plotSnapshots,
                  bool                 a_doLinINS,
                  bool                 a_doAdjoint,
                  bool                 a_doFirstOrderFreDeriv,
                  const EBIndexSpace* const a_ebisPtr)
{
  m_isDefined = true;
  m_params = a_params;
  m_baseflowIBCFact = a_baseflowIBC;
  m_solverIBCFact = a_solverIBC;
  m_coarsestDomain = a_coarsestDomain;
  m_viscosity = a_viscosity;
  m_ebisPtr = a_ebisPtr;
  m_refRatio = m_params.m_refRatio;
  m_coarsestDx = m_params.m_domainLength/Real(a_coarsestDomain.size(0));
  m_plotSnapshots = a_plotSnapshots;
  m_doLinINS = a_doLinINS;
  m_doAdjoint = a_doAdjoint;
  m_doFirstOrderFreDeriv = a_doFirstOrderFreDeriv;
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
//  int retval = SpaceDim + 1; // SpaceDim components of velocity + Pressure
  int retval = SpaceDim; // Only velocity comps; no pressure
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
/*
  { // copy pressure
    int ncomp = 1;
    Interval srcInterv(0,0);
    Interval dstInterv(SpaceDim, SpaceDim);
    LevelData<EBCellFAB> tempLD(a_levDBL, ncomp, IntVect::Unit*getnGhost(), ebcellfact);
    read<EBCellFAB>(a_handleIn, tempLD, "pres", a_levDBL, Interval(), false);
    tempLD.copyTo(srcInterv, *a_levBaseflow, dstInterv);
  }
*/
}
/*********/
#endif
/*********/
/*********/
bool EBAMRINSInterface::
isLinearSolver() const
{
  return m_doLinINS;
}
/*********/
void EBAMRINSInterface::
computeSolution(Epetra_Vector& a_y, const Epetra_Vector& a_x, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, const std::string& a_baseflowFile, double a_eps, double a_integrationTime, bool a_incOverlapData) const
{
  CH_assert(m_isDefined);

  int nlevels = a_baseflowDBL.size();
  int nComp = this->nComp();

  s_callCounter++;

  pout() << "doing y = EBAMRINSOp*x" << endl;

  if (m_doLinINS)
  {
    // make a_y = f(Ubar + eps*Uprime);
    // the solution should be independent of eps
    EBAMRLinINS solver(m_params, *m_baseflowIBCFact, *m_solverIBCFact, m_coarsestDomain, m_viscosity, m_doAdjoint);

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
//      const Vector<LevelData<EBCellFAB>* > presSoln = solver.getPresNew();

      int nVeloComp = veloSoln[0]->nComp();
//      int nPresComp = presSoln[0]->nComp();
      int totComp = this->nComp();
//    CH_assert(totComp == nVeloComp + nPresComp);

//      ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, veloSoln, 0., 1., 0, 0, nVeloComp, totComp, a_incOverlapData, m_refRatio);

      ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, veloSoln, 0., 1./a_eps, 0, 0, nVeloComp, totComp, a_incOverlapData, m_refRatio);

//    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, presSoln, 0., 1., nVeloComp, 0, nPresComp, totComp, a_incOverlapData, m_refRatio);
  } 

  else
  {
    {
      // compute Frechet Derivative 
  
      // make a_yStar = (f(Ubar + eps*Uprime))
      EBAMRNoSubcycle solver(m_params, *m_baseflowIBCFact, m_coarsestDomain, m_viscosity);
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
//      const Vector<LevelData<EBCellFAB>* > presSoln = solver.getPresNew();

      int nVeloComp = veloSoln[0]->nComp();
//      int nPresComp = presSoln[0]->nComp();
      int totComp = this->nComp();
//    CH_assert(totComp == nVeloComp + nPresComp);

      ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, veloSoln, 0., 1., 0, 0, nVeloComp, totComp, a_incOverlapData, m_refRatio);
//    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, presSoln, 0., 1., nVeloComp, 0, nPresComp, totComp, a_incOverlapData, m_refRatio);

    }

    if (m_doFirstOrderFreDeriv) // make a_y = a_yStar - f(Ubar)
    {
      Epetra_Vector baseflowVec(a_x);
      getBaseflow(baseflowVec, a_baseflowDBL, a_baseflowEBLG, a_baseflowFile, a_incOverlapData); 
      a_y.Update(-1., baseflowVec, 1.);
    }
    else // make a_y = a_yStar - (f(Ubar - eps*Uprime))
    {
      EBAMRNoSubcycle solver(m_params, *m_baseflowIBCFact, m_coarsestDomain, m_viscosity);
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
//      const Vector<LevelData<EBCellFAB>* > presSoln = solver.getPresNew();

      int nVeloComp = veloSoln[0]->nComp();
//      int nPresComp = presSoln[0]->nComp();
      int totComp = this->nComp();
//    CH_assert(totComp == nVeloComp + nPresComp);

      ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, veloSoln, 1., -1., 0, 0, nVeloComp, totComp, a_incOverlapData, m_refRatio);
//    ChomboEpetraOps::addChomboDataToEpetraVec(&a_y, presSoln, 1., -1., nVeloComp, 0, nPresComp, totComp, a_incOverlapData, m_refRatio);

    }

    double factor = (a_eps < 1.e-12) ? 1. : a_eps;
    double scale = m_doFirstOrderFreDeriv ? 1./factor : 0.5/factor;

    int checkScale = a_y.Scale(scale);
    CH_assert(checkScale == 0);

  } // end Frechet derivative

  pout() << "computed y = EMAMRINSOp*x" << endl;
  pout() << "Matrix-Vector Multiplication was computed " << s_callCounter << " times" << endl;

  if (m_plotSnapshots)
  {  
    std::string pltNameY = "vector_y_at_step"+SSTR(s_callCounter) + ".hdf5";
    std::string pltNameX = "vector_x_at_step"+SSTR(s_callCounter) + ".hdf5";

    pout() << "plotting a_x" << endl;
    plotEpetraVector(a_x, a_baseflowDBL, a_baseflowEBLG, pltNameX, a_incOverlapData,false);

    pout() << "plotting a_y" << endl;
    plotEpetraVector(a_y, a_baseflowDBL, a_baseflowEBLG, pltNameY, a_incOverlapData,false);
  }
}
/*********/
/*********/
void EBAMRINSInterface::
plotEpetraVector(const Epetra_Vector& a_v, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_fileName, bool a_incOverlapData, bool a_writeCheckpoint) const
{
  pout() << "plotting EigenEvector" << endl;

  if (!m_doLinINS)
  {
    EBAMRNoSubcycle solver(m_params, *m_baseflowIBCFact, m_coarsestDomain, m_viscosity);

    std::string baseflowFile = "dummyFile";

    solver.setupForStabilityRun(a_v, a_baseflowDBL, a_baseflowEBLG, baseflowFile, 1.0, a_incOverlapData, true);
//    solver.run(0,0);
    solver.concludeStabilityRun(&a_fileName, a_writeCheckpoint);
  }
  else
  {
    EBAMRLinINS solver(m_params, *m_baseflowIBCFact, *m_solverIBCFact, m_coarsestDomain, m_viscosity, m_doAdjoint);

    std::string baseflowFile = "dummyFile";

    solver.setupForStabilityRun(a_v, a_baseflowDBL, a_baseflowEBLG, baseflowFile, 1.0, a_incOverlapData, true);
//    solver.run(0,0);
    solver.concludeStabilityRun(&a_fileName, a_writeCheckpoint);
  }
}
/*********/
/*********/
void EBAMRINSInterface::
getBaseflow(Epetra_Vector& a_v, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_baseflowFile, bool a_incOverlapData) const
{
  Vector<LevelData<EBCellFAB>* > baseflow;
  int nlevels = a_baseflowDBL.size();
  baseflow.resize(nlevels);

  int nComp = this->nComp();
  int nGhost = getnGhost();

  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    EBCellFactory ebcellfact(a_baseflowEBLG[ilev].getEBISL());
    baseflow[ilev] = new LevelData<EBCellFAB>(a_baseflowDBL[ilev], nComp, nGhost*IntVect::Unit, ebcellfact);
  }

  HDF5Handle handleIn(a_baseflowFile, HDF5Handle::OPEN_RDONLY);
  readFileAndCopyToBaseflow(baseflow, a_baseflowFile, handleIn);
  handleIn.close();

  ChomboEpetraOps::rollChomboDataToEpetraVec(baseflow, &a_v, a_incOverlapData, getRefRatio());

  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    delete baseflow[ilev];
  }
}
/*********/
/*********/
void EBAMRINSInterface::
getBaseflow(Vector<LevelData<EBCellFAB>* >& a_baseLD, const Vector<DisjointBoxLayout>& a_baseflowDBL, const Vector<EBLevelGrid>& a_baseflowEBLG, std::string a_baseflowFile, bool a_incOverlapData) const
{
  HDF5Handle handleIn(a_baseflowFile, HDF5Handle::OPEN_RDONLY);
  readFileAndCopyToBaseflow(a_baseLD, a_baseflowFile, handleIn);
  handleIn.close();
}
/*********/
/*********/
