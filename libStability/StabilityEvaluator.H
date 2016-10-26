/*
 *
 *
 *
 *
 */

// This uses Trilinos package to compute eigen values

#ifndef _StabilityEvaluator_H_
#define _StabilityEvaluator_H_

#include "TrilinosSolverInterfaceFactory.H"
#include "TrilinosSolverInterface.H"

//includes from Trilinos
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Version.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>

using namespace std;

template <typename T>
class StabilityEvaluator
{

public:

  // setup stabilty evaluator:
  // a_solverInterfaceFactory creates the solverInterface to the solver
  // a_subspaceDim is the Krylov subspace dim
  // a_eps is the magnitude of the perturbation
  // a_timestep is the time between the snapshots that span Krylov subspace
  // a_baseFlowFile is the file containing the baseFlow

  // this constructs serial version (sets Epetra_Comm = Epetra_SerialComm)

  StabilityEvaluator(int                        a_subspaceDim, 
                     double                     a_eps, 
                     double                     a_timeStep,
                     string                     a_baseflowFile,
                     MPI_Comm*                  a_commPtr)
  {
    define(a_subspaceDim, a_eps, a_timeStep, a_baseflowFile, a_commPtr);
  }

  // setup stabilty evaluator:
  // a_solverInterfaceFactory creates the solverInterface to the solver
  // a_subspaceDim is the Krylov subspace dim
  // a_eps is the magnitude of the perturbation
  // a_timestep is the time between the snapshots that span Krylov subspace
  // a_baseFlowFile is the file containing the baseFlow
  // a_comm is the communication between processors: serial or parallel comm
  // if a_comm is NULL, it sets up for serial run else it uses solver's MPI_Comm
  void define(int       a_subspaceDim, 
              double    a_eps, 
              double    a_timeStep,
              string    a_baseflowFile,
              MPI_Comm* a_commPtr)
  {
    m_eps             = a_eps;
    m_subspaceDim     = a_subspaceDim;
    m_timeStep        = a_timeStep;
    m_baseflowFile    = a_baseflowFile;

    if (a_commPtr == NULL)
      {
        Epetra_SerialComm* commPtr = new Epetra_SerialComm ();
        m_commPtr = (Epetra_Comm*)commPtr;
        
        m_isParallelRun = false;
      }
    else
      {
        Epetra_MpiComm* commPtr = new Epetra_MpiComm (*a_commPtr);
        m_commPtr = (Epetra_Comm*)commPtr;
        m_isParallelRun = true;
      }

    m_isDefined       = true;
  }

  ~StabilityEvaluator()
  {
    delete m_commPtr;
    delete m_interface;
  }

  void
  exampleRoutine ()
  {
    if (m_commPtr->MyPID () == 0) {
      // On (MPI) Process 0, print out the Epetra software version.
      std::cout << Epetra_Version () << std::endl << std::endl;
    }
  }

private:

  bool m_isDefined;
  bool m_isParallelRun;

  int     m_subspaceDim;
  double  m_eps;
  double  m_timeStep;
  string  m_baseflowFile;  

  Epetra_Comm* m_commPtr;

  TrilinosSolverInterface<T>* m_interface;

  // disallowed for all usual reasons
  StabilityEvaluator()
  {;} 
};

#endif
