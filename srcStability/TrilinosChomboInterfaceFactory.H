/*
 *
 *
 *
 *
 */

#ifndef _Trilinos_Chombo_Interface_Factory_H_
#define _Trilinos_Chombo_Interface_Factory_H_

#include "TrilinosSolverInterfaceFactory.H"
#include "TrilinosChomboInterface.H"
#include "MayDay.H"

#include <Teuchos_RCPDecl.hpp>

class TrilinosChomboInterfaceFactory : public TrilinosSolverInterfaceFactory
{
public:

  TrilinosChomboInterfaceFactory()
  {
     MayDay::Error("TrilinosChomboInterfaceFactory::invalid constructor");
  }

  TrilinosChomboInterfaceFactory(RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact)
  {
    m_solverInterfaceFact = a_solverInterfaceFact;
  }

  virtual ~TrilinosChomboInterfaceFactory()
  {;}

  Teuchos::RCP<TrilinosSolverInterface> create() const
  {
    Teuchos::RCP<TrilinosChomboInterface> ChomboInterface = Teuchos::rcp (new TrilinosChomboInterface (m_solverInterfaceFact) );

    Teuchos::RCP<TrilinosSolverInterface> retval = static_cast<Teuchos::RCP<TrilinosSolverInterface> > (ChomboInterface);

    return retval;
  }

protected:

private:

  RefCountedPtr<ChomboSolverInterfaceFactory> m_solverInterfaceFact;
  
};

#endif
