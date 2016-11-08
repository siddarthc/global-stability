/*
 *
 *
 *
 *
 */

#ifndef _CHOMBO_SOLVER_INTERFACE_FACTORY_H_
#define _CHOMBO_SOLVER_INTERFACE_FACTORY_H_

#include "ChomboSolverInterface.H"

class ChomboSolverInterfaceFactory
{
public:

  ChomboSolverInterfaceFactory()
  {;}

  virtual ~ChomboSolverInterfaceFactory()
  {;}

  virtual ChomboSolverInterface* create() const = 0; 

};

#endif
