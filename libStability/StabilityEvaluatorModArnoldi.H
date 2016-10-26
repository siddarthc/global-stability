/*
 *
 *
 *
 *
 */

// implementation of the Modified Arnoldi algorithm
// Incomplete!!!

#ifndef _StabilityEvaluator_H_
#define _StabilityEvaluator_H_

#include "SolverInterfaceFactory.H"
#include "SolverInterface.H"

using namespace std;

template <typename T>
class StabilityEvaluator
{

public:

/*
  StabilityEvaluator()
  {

  }
*/

  // setup stabilty evaluator:
  // a_solverInterfaceFactory creates the solverInterface to the solver
  // a_subspaceDim is the Krylov subspace dim
  // a_eps is the magnitude of the perturbation
  // a_timestep is the time between the snapshots that span Krylov subspace
  // a_baseFlowFile is the file containing the baseFlow

  StabilityEvaluator(SolverInterfaceFactory<T>& a_solverInterfaceFact, 
                     int                        a_subspaceDim, 
                     double                     a_eps, 
                     double                     a_timeStep,
                     string                     a_baseflowFile)
  {
    define(a_solverInterfaceFact, a_subspaceDim, a_eps, a_timeStep, a_baseflowFile);
  }

  // setup stabilty evaluator:
  // a_solverInterfaceFactory creates the solverInterface to the solver
  // a_subspaceDim is the Krylov subspace dim
  // a_eps is the magnitude of the perturbation
  // a_timestep is the time between the snapshots that span Krylov subspace
  // a_baseFlowFile is the file containing the baseFlow
  
  void define(SolverInterfaceFactory<T>& a_solverInterfaceFact, 
              int                        a_subspaceDim, 
              double                     a_eps, 
              double                     a_timeStep,
              string                     a_baseflowFile)
  {
    m_solverInterface = a_solverInterfaceFact.create();
    m_eps             = a_eps;
    m_subspaceDim     = a_subspaceDim;
    m_timeStep        = a_timeStep;
    m_baseflowFile    = a_baseflowFile;
    m_isDefined       = true;

    allocate();
  }

  ~StabilityEvaluator()
  {
    for(int i = 0; i < m_subspaceDim; i++)
      {
        delete [] m_hessenMatrix[i];
      }

    delete [] m_hessenMatrix;
  }


  void computeStability()
  {
    m_solverInterface->setupForStabilityRun(m_baseflowFile);
    m_solverInterface->copySolution(m_baseflow);
    fillInitialKrylovVector();
    computeArnoldiIterations();
    clearKrylovVectors();
       
//    computeEigenModes();

  }

private:

  void fillInitialVector()
  {
    m_solverInterface->fillWithNoise(m_krylovVector[0]);
    m_solverInterface->normalize(m_krylovVector[0]);
  }
  
  void allocate()
  {
    m_hessenMatrix  = new double*[m_subspaceDim];
    for (int i = 0; i < m_subspaceDim; i++)
      {
        m_hessenMatrix[i] = new double[m_subspaceDim];
      }

    m_krylovVectors = new T[m_subspaceDim];
  }

  void computeArnoldiIterations()
  {
    // Arnoldi iteration loop
    for (int i = 0; i < m_subspaceDim; i++)
      {
        m_solverInterface->applyLinearOp(m_baseflow, m_krylovVectors[i], m_eps, m_timeStep);
      }
  }

  void clearKrylovVectors()
  {
    delete [] m_krylovVectors;
  }

  bool m_isDefined;

  SolverInterface<T>* m_solverInterface;

  int     m_subspaceDim;
  double  m_eps;
  double  m_timeStep;
  string  m_baseflowFile;  

  double** m_hessenMatrix; // a 2D array
  T*       m_krylovVectors;
  T        m_baseflow;  

  // disallowed for all usual reasons
  StabilityEvaluator()
  {;} 
};

#endif
