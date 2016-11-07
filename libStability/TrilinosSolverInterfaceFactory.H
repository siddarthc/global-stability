/*
 *
 *
 *
 *
 */

#ifndef _TrilinosSolverInterfaceFactory_H_
#define _TrilinosSolverInterfaceFactory_H_

#include <Teuchos_RCPDecl.hpp>
#include "TrilinosSolverInterface.H"


class TrilinosSolverInterfaceFactory
{
public:

  TrilinosSolverInterfaceFactory()
  {

  }

  virtual ~TrilinosSolverInterfaceFactory()
  {

  }

  virtual Teuchos::RCP<TrilinosSolverInterface> create() const = 0;
};

#endif
