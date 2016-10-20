/*
 *
 *
 *
 *
 */

#ifndef _SolverInterfaceFactory_H_
#define _SolverInterfaceFactory_H_

#include "SolverInterface.H"

template <typename T>
class SolverInterfaceFactory
{
public:

  SolverInterfaceFactory()
  {;}

  virtual ~SolverInterfaceFactory()
  {;}

  virtual SolverInterface<T>* create() const = 0;
};

#endif
