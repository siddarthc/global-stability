/*
 *
 *
 *
 */

#ifndef _ANASAZI_OPERATOR_SOLVER_EXT_H_
#define _ANASAZI_OPERATOR_SOLVER_EXT_H_

#include <Epetra_MultiVector.h>
#include <AnasaziOperator.hpp>
#include "TrilinosSolverInterface.H"
#include <Teuchos_RCPDecl.hpp>

using Teuchos::RCP;

namespace Anasazi {

  class OperatorSolverExt : public Operator<double>
  {
  public:

    OperatorSolverExt(const RCP<TrilinosSolverInterface>& a_solverInterface);

    void Apply(const MultiVec<double>& x, MultiVec<double>& y) const;

  private:
    RCP<TrilinosSolverInterface> m_solverInterface;
  };

/*
  template<>
  class OperatorTraits < double, MultiVec<double>, OperatorSolverExt >
  {
  public:

    static void Apply(const AnasaziOperatorSolverExt& Op,
                      const Epetra_MultiVector& x,
                      Epetra_MultiVector& y)
    {
      Op.Apply( x, y );
    }

  };
*/
}

#endif
