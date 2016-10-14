/*
 *
 *
 *
 *
 *
 */

#ifndef _SFDOpFunctions_H_
#define _SFDOpFunctions_H_


template <typename T>
class SFDOpFunctions
{

public:

  SFDOpFunctions()
  {

  }

  ~SFDOpFunctions()

  {

  }

  virtual void aXbY(T&            a_result,
                    const T&      a_X,
                    const double& a_aCoef,
                    const T&      a_Y,
                    const double& a_bCoef) = 0;

  virtual void copierFunc(T&       a_clone,
                          const T& a_data) = 0;
};

#endif
