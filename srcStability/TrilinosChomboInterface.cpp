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
TrilinosChomboInterface(const RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact, bool a_incOverlapData, bool a_doWeighting)
{
  m_solverInterface   = a_solverInterfaceFact->create();
  m_isBaseflowSet     = false;
  m_isVolWeightsSet   = false;
  m_domainVolume      = 0.;
  m_volWeights        = NULL;
  m_incOverlapData    = a_incOverlapData;
  m_doWeighting       = a_doWeighting;
  m_isBaseflowNormSet = false;
  m_baseflowL2Norm    = -1.;
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
setBaseflow(const std::string& a_baseflowFile, const Epetra_Comm* a_commPtr)
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

  // volWeights:
  if (m_doWeighting)
  {
    Epetra_Map map = getEpetraMap(a_commPtr);
    m_volWeights = new Epetra_Vector(map);
    ChomboEpetraOps::getVolWeights(m_volWeights, m_domainVolume, m_baseflowDBL, m_baseflowEBLG, m_solverInterface->nComp(), m_solverInterface->getCoarsestDx(), m_solverInterface->getRefRatio(), m_incOverlapData);
    m_isVolWeightsSet = true;
  }

  computeBaseflowNorm(a_commPtr);
}

/*********/
int TrilinosChomboInterface::
nElementsOnThisProc() const
{
  CH_assert(isSetupForStabilityRun());
  int nComp = m_solverInterface->nComp();
  int retval = ChomboEpetraOps::getnElementsOnThisProc(m_baseflowDBL, m_baseflowEBLG, nComp, m_incOverlapData, m_solverInterface->getRefRatio());
  return retval;
}

/*********/
Epetra_Map TrilinosChomboInterface::
getEpetraMap(const Epetra_Comm* a_commPtr) const
{
  CH_assert(isSetupForStabilityRun());
  int nComp = m_solverInterface->nComp();
  Epetra_Map retval = ChomboEpetraOps::getEpetraMap(m_baseflowDBL, m_baseflowEBLG, nComp, a_commPtr, m_incOverlapData, m_solverInterface->getRefRatio());
  return retval;
}

/*********/
int TrilinosChomboInterface::
computeL2Norm(const Epetra_Vector& a_v, double& a_result) const
{
  CH_assert(isSetupForStabilityRun());
  int retval;
  if (!m_doWeighting)
  {
    retval = ChomboEpetraOps::computeL2Norm(a_result, a_v);
  }
  else
  {
    CH_assert(m_isVolWeightsSet);
    CH_assert(m_volWeights != NULL);
    retval = ChomboEpetraOps::computeWeightedL2Norm(a_result, a_v, *m_volWeights);
  }

  return retval;
}

/*********/
int TrilinosChomboInterface::
computeDotProd(const Epetra_Vector& a_v1, const Epetra_Vector& a_v2, double& a_result) const
{
  CH_assert(isSetupForStabilityRun());
  int retval;
  if (!m_doWeighting)
  {
    retval = ChomboEpetraOps::computeDotProduct(a_result, a_v1, a_v2);
  }
  else
  {
    CH_assert(m_isVolWeightsSet);
    CH_assert(m_volWeights != NULL);
    retval = ChomboEpetraOps::computeWeightedDotProduct(a_result, a_v1, a_v2, *m_volWeights);
  }

  return retval;
}

/*********/
void TrilinosChomboInterface::
computeSolution(const Epetra_Vector& a_x, Epetra_Vector& a_y) const
{
  CH_assert(isSetupForStabilityRun());
  CH_assert(m_isBaseflowNormSet);
  CH_assert(m_baseflowL2Norm > 0.);

  double vecNorm;
  computeL2Norm(a_x, vecNorm);

  bool isLinearSolver = m_solverInterface->isLinearSolver();

  double pertSize;  

  if (isLinearSolver)
  {
    pertSize = m_eps;
  }
 
  else
  {
    pertSize = m_eps*m_baseflowL2Norm/vecNorm; // Theofillis paper
//    double pertSize = sqrt((1 + m_baseflowL2Norm)*m_eps)/vecNorm; // refer: JFNK: survey of applications... by Knoll and Keyes
//    double pertSize = m_eps*m_baseflowL2Norm; // refer: https://arxiv.org/pdf/1502.03701.pdf
  }

  pout() << "Base flow L2 norm = " << m_baseflowL2Norm << endl;
  pout() << "vector L2 norm = " << vecNorm << endl;
  pout() << "pert size = " << pertSize << endl;
  pout() << "domain volume = " << m_domainVolume << endl;

  m_solverInterface->computeSolution(a_y, a_x, m_baseflowDBL, m_baseflowEBLG, m_baseflowFile, pertSize, m_integrationTime, m_incOverlapData);
}
/*********/
void TrilinosChomboInterface::
plotEpetraVector(const Epetra_Vector& a_v, std::string a_plotName) const
{
  CH_assert(isSetupForStabilityRun());
  std::string name = a_plotName + ".hdf5";
  m_solverInterface->plotEpetraVector(a_v, m_baseflowDBL, m_baseflowEBLG, name, m_incOverlapData);
}
/*********/
void TrilinosChomboInterface::
computeBaseflowNorm(const Epetra_Comm* a_commPtr)
{
  CH_assert(m_isBaseflowSet)
  
  Epetra_Map tmpMap = getEpetraMap(a_commPtr);
  Epetra_Vector tmpVector(tmpMap);

  m_solverInterface->getBaseflow(tmpVector, m_baseflowDBL, m_baseflowEBLG, m_baseflowFile, m_incOverlapData);

  computeL2Norm(tmpVector, m_baseflowL2Norm); 
  m_isBaseflowNormSet = true;
}

/*********/
