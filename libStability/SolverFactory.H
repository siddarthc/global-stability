/*
 *
 *
 *
 *
 */

#ifndef _SolverFactory_H_
#define _SolverFactory_H_

#include "Solver.H"

class SolverFactory
{
public:

  SolverFactory()
  {;}

  virtual ~SolverFactory()
  {;}

  virtual Solver* create() const = 0;
};

#endif
