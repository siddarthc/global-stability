/*
 *
 *
 *
 *
 */

#include "TrilinosChomboInterface.H"
#include "ChomboEpetraOps.H"
#include "EBEllipticLoadBalance.H"
#include "EBLevelGrid.H"
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

  int nlevels = m_baseflow.size();
  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    delete m_baseflow[ilev];
  }
}

/*********/
void TrilinosChomboInterface::
setBaseflow(std::string a_baseflowFile)
{

#ifdef CH_USE_HDF5

  HDF5Handle handleIn(a_baseflowFile, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handleIn);

  int finestLevel = header.m_int["finest_level"];

  int finestLevelFromParmParse = m_solverInterface->getMaxLevelFromParmParse();

  if (finestLevel != finestLevelFromParmParse)
  {
     MayDay::Warning("TrilinosChomboInterface::setBaseflow - finestLevel from baseflow file is inconsistent with finestLevel from input file. Proceeding with finestLevel from baseflow file");
  }

  m_baseflow.resize(finestLevel+1, NULL);

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
    DisjointBoxLayout levelGrids = DisjointBoxLayout(vboxGrids,proc_map);

    EBLevelGrid levelEBLG = EBLevelGrid(levelGrids, levelDomain, nEBGhost, ebisPtr);
    EBISLayout levelEBISL = levelEBLG.getEBISL();
 
    // define baseflow
    EBCellFactory ebcellfact(levelEBISL);
    m_baseflow[ilev] = new LevelData<EBCellFAB>(levelGrids, nComp, nGhost*IntVect::Unit, ebcellfact);
    EBLevelDataOps::setToZero(*(m_baseflow[ilev]));
    m_solverInterface->readFileAndCopyToBaseflow(m_baseflow[ilev], levelGrids, levelEBISL, levelDomain, a_baseflowFile, handleIn);
  }

//  EBAMRDataOps::setToZero(m_baseflow);
//  m_solverInterface->readFileAndCopyToBaseflow(m_baseflow, a_baseflowFile, handleIn);
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
  int retval = ChomboEpetraOps::getnElementsOnThisProc(m_baseflow);
  return retval;
}

/*********/
Epetra_Map TrilinosChomboInterface::
getEpetraMap(const Epetra_Comm* a_commPtr) const
{
  CH_assert(isSetupForStabilityRun());
  Epetra_Map retval = ChomboEpetraOps::getEpetraMap(m_baseflow, a_commPtr);
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
  m_solverInterface->computeSolution(a_y, a_x, m_baseflow, m_eps, m_integrationTime);
}
