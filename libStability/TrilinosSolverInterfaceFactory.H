/*
 *
 *
 *
 *
 */

#ifndef _TrilinosSolverInterfaceFactory_H_
#define _TrilinosSolverInterfaceFactory_H_

#include "TrilinosSolverInterface.H"

template <typename T>
class TrilinosSolverInterfaceFactory
{
public:

  TrilinosSolverInterfaceFactory()
  {;}

  virtual ~TrilinosSolverInterfaceFactory()
  {;}

  virtual TrilinosSolverInterface<T>* create() const = 0;
};

#endif
