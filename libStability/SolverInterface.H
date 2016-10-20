/*
 *
 *
 *
 *
 */

#ifndef _SolverInterface_H_
#define _SolverInterface_H_

#include <string>

template <typename T>
class SolverInterface
{
public:

  SolverInterface()
  {

  }

  virtual ~SolverInterface()
  {

  }

  virtual void setupForStabilityRun(std::string a_string) = 0;

  virtual void copySolution(T& a_state) = 0;

  virtual void fillWithNoise(T& a_state) = 0;

  virtual void normalize(T& a_state) = 0;

  virtual void run(double a_startTime, double a_endTime) = 0;
};

#endif
