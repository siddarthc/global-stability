/*
 *
 *
 *
 *
 */

#include "TrilinosChomboInterface.H"

/*********/
TrilinosChomboInterface::
TrilinosChomboInterface(const RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact)
{
  m_solverInterface = a_solverInterfaceFact->create();
}

/*********/
TrilinosChomboInterface::
~TrilinosChomboInterface()
{
  delete m_solverInterface;
}

/*********/
void TrilinosChomboInterface::
setupForStabilityRun(std::string a_baseflowFile)
{

}

/*********/
int TrilinosChomboInterface::
nElementsOnThisProc() const
{
  return 0;
}

/*********/
int TrilinosChomboInterface::
computeL2Norm(const Epetra_MultiVector& a_mv, double* a_result) const
{
  return 1;
}

/*********/
int TrilinosChomboInterface::
computeDotProd(const Epetra_MultiVector& a_mv1, const Epetra_MultiVector& a_mv2, double* a_result) const
{
  return 1;
}

/*********/
void TrilinosChomboInterface::
computeSolution(const Epetra_MultiVector& a_x, Epetra_MultiVector& a_y) const
{

}
