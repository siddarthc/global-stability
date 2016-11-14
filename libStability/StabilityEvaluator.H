/*
 *
 *
 *
 *
 */

// This uses Trilinos package to compute eigenvalues

#ifndef _StabilityEvaluator_H_
#define _StabilityEvaluator_H_

#include "TrilinosSolverInterfaceFactory.H"

//includes from Trilinos
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_RCPDecl.hpp>
// Include header for block Krylov-Schur solver
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
// Include header to define basic eigenproblem Ax = \lambda*Bx
#include "AnasaziBasicEigenproblem.hpp"
// Include header to provide Anasazi with Epetra adapters.  If you
// plan to use Tpetra objects instead of Epetra objects, include
// AnasaziTpetraAdapter.hpp instead; do analogously if you plan to use
// Thyra objects instead of Epetra objects.

#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziEpetraSolverAdapter.H"
#include "AnasaziOperatorSolverExt.H"

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;

//typedef Epetra_MultiVector MV;
//typedef Anasazi::EpetraMultiVec MV;
typedef Teuchos::ScalarTraits<double>                   SCT;
typedef SCT::magnitudeType                              MagnitudeType;
typedef Anasazi::EpetraMultiVecSolverExt                MV;
typedef Anasazi::OperatorSolverExt                      OP;
typedef Anasazi::MultiVec<double>                       BaseMV;
typedef Anasazi::Operator<double>                       BaseOP;
typedef Anasazi::MultiVecTraits<double, BaseMV >        MVT;
typedef Anasazi::OperatorTraits<double, BaseMV, BaseOP> OPT;

class StabilityEvaluator
{

public:

  // setup stabilty evaluator:
  // a_solverInterfaceFactory creates the solverInterface to the solver
  // a_subspaceDim is the Krylov subspace dim
  // a_eps is the magnitude of the perturbation
  // a_timestep is the time between the snapshots that span Krylov subspace
  // a_baseflowFile is the name of the file which contains the baseflow
  // a_solverInterfaceFact is the factory that produces solverInterface when needed
  // if a_comm is NULL, it sets up for serial run else it uses solver's MPI_Comm

  StabilityEvaluator(int                                 a_subspaceDim, 
                     double                              a_eps, 
                     double                              a_timeStep,
                     string                              a_baseflowFile,
                     RCP<TrilinosSolverInterfaceFactory> a_solverInterfaceFact,
                     MPI_Comm*                           a_commPtr);

  // setup stabilty evaluator:
  // a_solverInterfaceFactory creates the solverInterface to the solver
  // a_subspaceDim is the Krylov subspace dim
  // a_eps is the magnitude of the perturbation
  // a_timestep is the time between the snapshots that span Krylov subspace
  // a_baseflowFile is the name of file which contains the baseflow
  // a_solverInterfaceFact is the factory that produces solverInterface when needed
  // a_comm is the communication between processors: serial or parallel comm
  // if a_comm is NULL, it sets up for serial run else it uses solver's MPI_Comm

  void define(int                                 a_subspaceDim, 
              double                              a_eps, 
              double                              a_timeStep,
              string                              a_baseflowFile,
              RCP<TrilinosSolverInterfaceFactory> a_solverInterfaceFact,
              MPI_Comm*                           a_commPtr);

  virtual ~StabilityEvaluator();

  // a_tol is tolerence for the eigen value solver
  // a_nev is number of eigen values to be computed
  // a_blockSize is the size of the block of vectors
  // a_maxRestarts is the max number of restarts allowed
  // a_which is tells how to sort the eigen values. Options:
  //   SM -> Sort in increasing order of magnitude
  //   SR -> Sort in increasing order of real part
  //   SI -> Sort in increasing order of imaginary part
  //   LM -> Sort in decreasing order of magnitude
  //   LR -> Sort in decreasing order of real part
  //   LI -> Sort in decreasing order of imaginary part
  // a_isOpSymmetric is true if Operator is symmetric
  // a_verbose: if true, prints eigen values

  // returns 0 if success

  int computeDominantModes(double a_tol,
                           int    a_nev,
                           int    a_numBlocks,
                           int    a_blockSize,
                           int    a_maxRestarts,
                           string a_which,
                           bool   a_isOpSymmetric,
                           bool   a_verbose);

  void
  exampleRoutine ();

private:

  bool m_isDefined;
  bool m_isParallelRun;

  int    m_subspaceDim;
  double m_eps;
  double m_integrationTime;
  string m_baseflowFile; 

  Epetra_Comm* m_commPtr;

  RCP<TrilinosSolverInterface> m_solverInterface;

  // disallowed for all usual reasons
  StabilityEvaluator()
  {;} 
};

#endif
