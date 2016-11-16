/*
 *
 *
 *
 *
 */

#include "TrilinosChomboInterface.H"
#include "ChomboEpetraOps.H"
#include "EBEllipticLoadBalance.H"
#include "EBCellFactory.H"
#include "EBAMRDataOps.H"
#include "CH_HDF5.H"

/*********/
TrilinosChomboInterface::
TrilinosChomboInterface(const RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact)
{
  m_solverInterface = a_solverInterfaceFact->create();
  m_isBaseflowSet   = false;
  m_isVolWeightsSet = false;
  m_domainVolume    = 0.;
  m_volWeights      = NULL;
}

/*********/
TrilinosChomboInterface::
~TrilinosChomboInterface()
{
  delete m_solverInterface;

  if (m_volWeights != NULL) delete m_volWeights;
}

/*********/
void TrilinosChomboInterface::
setBaseflow(const std::string& a_baseflowFile)
{
  m_baseflowFile = a_baseflowFile;

#ifdef CH_USE_HDF5

  HDF5Handle handleIn(m_baseflowFile, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handleIn);

  int finestLevel = header.m_int["finest_level"];

  int finestLevelFromParmParse = m_solverInterface->getMaxLevelFromParmParse();

  if (finestLevel != finestLevelFromParmParse)
  {
     MayDay::Warning("TrilinosChomboInterface::setBaseflow - finestLevel from baseflow file is inconsistent with finestLevel from input file. Proceeding with finestLevel from baseflow file");
  }

//  m_baseflow.resize(finestLevel+1, NULL);

  m_baseflowDBL.resize(finestLevel+1);
  m_baseflowEBLG.resize(finestLevel+1);

  int nGhost = m_solverInterface->getnGhost();
  int nEBGhost = m_solverInterface->getnEBGhost();
  const EBIndexSpace* ebisPtr = m_solverInterface->getEBISPtr();
  int nComp = m_solverInterface->nComp();

  // get all the grids
  for (int ilev=0; ilev <= finestLevel; ilev++)
  {
    handleIn.setGroupToLevel(ilev);
    // Get the grids
    Vector<Box> vboxGrids;
    const int gridStatus = read(handleIn, vboxGrids);
    if (gridStatus != 0)
    {
      MayDay::Error("readCheckpointLevel: file has no grids");
    }

    Vector<int> proc_map;
    ProblemDomain levelDomain;
    m_solverInterface->getLevelDomain(&levelDomain, ilev);    

    EBEllipticLoadBalance(proc_map, vboxGrids, levelDomain, false, ebisPtr);
    m_baseflowDBL[ilev] = DisjointBoxLayout(vboxGrids, proc_map);
    m_baseflowEBLG[ilev] = EBLevelGrid(m_baseflowDBL[ilev], levelDomain, nEBGhost, ebisPtr);

  }

  m_isBaseflowSet = true;

  handleIn.close();
#else

  MayDay::Error("Chombo needs HDF5 to read baseflowFile");

#endif

}

/*********/
int TrilinosChomboInterface::
nElementsOnThisProc() const
{
  CH_assert(isSetupForStabilityRun());
  int nComp = m_solverInterface->nComp();
  int retval = ChomboEpetraOps::getnElementsOnThisProc(m_baseflowDBL, m_baseflowEBLG, nComp);
  return retval;
}

/*********/
Epetra_Map TrilinosChomboInterface::
getEpetraMap(const Epetra_Comm* a_commPtr) const
{
  CH_assert(isSetupForStabilityRun());
  int nComp = m_solverInterface->nComp();
  Epetra_Map retval = ChomboEpetraOps::getEpetraMap(m_baseflowDBL, m_baseflowEBLG, nComp, a_commPtr);
  return retval;
}

/*********/
int TrilinosChomboInterface::
computeL2Norm(const Epetra_Vector& a_v, double& a_result) const
{
  CH_assert(isSetupForStabilityRun());
  int retval = ChomboEpetraOps::computeL2Norm(a_result, a_v);
  return retval;
}

/*********/
int TrilinosChomboInterface::
computeDotProd(const Epetra_Vector& a_v1, const Epetra_Vector& a_v2, double& a_result) const
{
  CH_assert(isSetupForStabilityRun());
  int retval = ChomboEpetraOps::computeDotProduct(a_result, a_v1, a_v2);
  return retval;
}

/*********/
void TrilinosChomboInterface::
computeSolution(const Epetra_Vector& a_x, Epetra_Vector& a_y) const
{
  CH_assert(isSetupForStabilityRun());
  m_solverInterface->computeSolution(a_y, a_x, m_baseflowDBL, m_baseflowEBLG, m_baseflowFile, m_eps, m_integrationTime);
}
/*********/
