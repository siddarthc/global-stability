/*
 *
 *
 *
 *
 */

#include "AnasaziOperatorSolverExt.H"
#include "AnasaziEpetraSolverAdapter.H"

namespace Anasazi {

  OperatorSolverExt::OperatorSolverExt(const RCP<TrilinosSolverInterface>& a_solverInterface)
  {
    m_solverInterface = a_solverInterface;
  }

  void OperatorSolverExt::Apply(const MultiVec<double>& x, MultiVec<double>& y) const
  {
    EpetraMultiVecSolverExt *x_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(x));
    TEUCHOS_TEST_FOR_EXCEPTION( x_vec==NULL,  std::invalid_argument, "Anasazi::OperatorSolverExt::Apply cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");

    EpetraMultiVecSolverExt *y_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&y);
    TEUCHOS_TEST_FOR_EXCEPTION( y_vec==NULL,  std::invalid_argument, "Anasazi::OperatorSolverExt::Apply cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");
    m_solverInterface->computeSolution(*x_vec, *y_vec);
  }
}
