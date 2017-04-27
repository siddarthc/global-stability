/*
 *
 *
 *
 *
 */

#ifndef _SFDOp_H_
#define _SFDOp_H_

#include <vector>
#include <cassert>
#include <cmath>
#include <memory>
#include <iostream>

#include "SFDOpFunctions.H"

using namespace std;

template <typename T>
class SFDOp
{
//  typedef shared_ptr<T> Tptr;
  typedef T* Tptr;

  typedef void (SFDOpFunctions<T>::*aXbYFunc_)(T&, const T&, const double&, const T&, const double&);
  typedef void (SFDOpFunctions<T>::*copierFunc_)(T&, const T&);

public:

  SFDOp() 
  {
    m_isDefined = false;
  }

// The damping scheme is defined by S(nFilters, smallestFilter, largestFilter): Refer Massa, JFM
// In addition, overall controlCoef can also be set as parameter
// if nFilters = 1, set smallest filter; value of largest filter doesn't matter
// if fixed grids, set extrapToNewGridFunc = NULL

  SFDOp(int                        a_nFilters, 
        double                     a_smallestFilter,
        double                     a_largestFilter,
        double                     a_controlCoef)
  {
    define(a_nFilters, a_smallestFilter, a_largestFilter, a_controlCoef);
  }

  void define(int                       a_nFilters,
              double                    a_smallestFilter,
              double                    a_largestFilter,
              double                    a_controlCoef)
  {
    //    m_counter             = 0;
    m_nFilters            = a_nFilters;

    setSFDParams(a_nFilters, a_smallestFilter, a_largestFilter, a_controlCoef);
  }

  ~SFDOp() 
  {
  }

  void operator() (T& a_q, vector<T*>& a_qBar, double a_dt, aXbYFunc_ a_aXbYFunc, copierFunc_ a_copier, SFDOpFunctions<T>* a_functions)
  {
    applyOperator(a_q, a_qBar, a_dt, a_aXbYFunc, a_copier, a_functions);
  }

  vector<double> getFilterWidth()
  {
    return m_filterWidth;
  }

  vector<double> getControlCoef()
  {
    return m_controlCoef;
  }

protected:

  void applyOperator(T& a_q, vector<T*>& a_qBar, double a_dt, aXbYFunc_ a_aXbYFunc, copierFunc_ a_copier, SFDOpFunctions<T>* a_functions)
  {
    for (int i = 0; i < m_nFilters; i++)
      {
        applyOperator(a_q, *a_qBar[i], i, a_dt, a_aXbYFunc, a_copier, a_functions);
      }
  }
  
  void applyOperator(T& a_q, T& a_qBar, int& a_index, double a_dt, aXbYFunc_ a_aXbYFunc, copierFunc_ a_copier, SFDOpFunctions<T>* a_functions)
  {
    T qOld, qBarOld;
    (a_functions->*a_copier)(qOld, a_q);
    (a_functions->*a_copier)(qBarOld, a_qBar);

    double X = a_dt*m_controlCoef[a_index]; 
    double Delta = m_filterWidth[a_index]/a_dt;
    double coeff = 1./(1. + X*Delta);

    double aCoef = coeff*(1. + X*Delta*exp(-(X + 1./Delta)));
    double bCoef = coeff*(X*Delta*(1. - exp(-(X + 1./Delta))));
    double cCoef = coeff*(1. - exp(-(X + 1./Delta)));
    double dCoef = coeff*(X*Delta + exp(-(X + 1./Delta)));

    (a_functions->*a_aXbYFunc)(a_q, qOld, aCoef, qBarOld, bCoef);
    (a_functions->*a_aXbYFunc)(a_qBar, qOld, cCoef, qBarOld, dCoef);
  }

  void setSFDParams(int a_nFilters, double a_smallestFilter, double a_largestFilter, double a_controlCoef)
  {
    m_isDefined   = true;
    m_filterWidth.resize(a_nFilters);
    m_controlCoef.resize(a_nFilters);

    if (a_nFilters == 1)
      {
        m_filterWidth[0] = a_smallestFilter;
        m_controlCoef[0] = a_controlCoef;
      }
    else
      {
        double ratio = a_largestFilter/a_smallestFilter;
        double denom = m_nFilters - 1;

        for (int i = 0; i < a_nFilters; i++)
          {
            m_filterWidth[i] = a_smallestFilter*pow(ratio, ((i)/denom));
            m_controlCoef[i] = a_controlCoef/a_nFilters;
          }
      }
  }

/*
  vector<double> getFilterWidth()
  {
    return m_filterWidth;
  }

  vector<double> getControlCoef()
  {
    return m_controlCoef;
  }
*/

protected:

  int m_nFilters;
  
  vector<double> m_filterWidth;
  vector<double> m_controlCoef;

private:

//  int  m_counter;
  bool m_isDefined;
  bool m_isInitialized;

};

#endif
