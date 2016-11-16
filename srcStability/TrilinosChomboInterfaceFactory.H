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

/*********/

  TrilinosChomboInterfaceFactory(RefCountedPtr<ChomboSolverInterfaceFactory>& a_solverInterfaceFact, bool a_incOverlapData)
  {
    m_solverInterfaceFact = a_solverInterfaceFact;
    m_incOverlapData      = a_incOverlapData;
  }

/*********/

  virtual ~TrilinosChomboInterfaceFactory()
  {;}

/*********/

  Teuchos::RCP<TrilinosSolverInterface> create() const
  {
    Teuchos::RCP<TrilinosChomboInterface> ChomboInterface = Teuchos::rcp (new TrilinosChomboInterface (m_solverInterfaceFact, m_incOverlapData) );

    Teuchos::RCP<TrilinosSolverInterface> retval = static_cast<Teuchos::RCP<TrilinosSolverInterface> > (ChomboInterface);

    return retval;
  }

/*********/

protected:

private:

  RefCountedPtr<ChomboSolverInterfaceFactory> m_solverInterfaceFact;
  
  bool m_incOverlapData;
};

#endif
