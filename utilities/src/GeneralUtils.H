/*
General Utility class - 19 April 2017

*/



#ifndef _GENERAL_UTILS_H_
#define _GENERAL_UTILS_H_

#include "ProblemDomain.H"
#include "EBIndexSpace.H"
#include "Vector.H"
#include "REAL.H"

class GeneralUtils
{
public:


  static void
  computeKappaDotProduct(Real& a_out,
                         const Vector<LevelData<EBCellFAB>* >& a_in1,
                         const Vector<LevelData<EBCellFAB>* >& a_in2,
                         Vector<int> a_refRatio,
                         Real a_coarsestDx);

  static void
  computeL2Norm(Real&                                 a_norm,
                const Vector<LevelData<EBCellFAB>* >& a_q,
	        const Vector<Real>&                   a_dxVect,
                const Vector<int>&                    a_refRatio,
                bool                                  a_incOverlapData = false);

  static void
  computeLinfNorm(Real&                                 a_norm,
                  const Vector<LevelData<EBCellFAB>* >& a_q,
                  const Vector<Real>&                   a_dxVect,
                  const Vector<int>&                    a_refRatio,
                  bool                                  a_incOverlapData = false);


  static void
  computeL0Norm(Real&                                 a_norm,
                  const Vector<LevelData<EBCellFAB>* >& a_q,
                  const Vector<Real>&                   a_dxVect,
                  const Vector<int>&                    a_refRatio,
                  bool                                  a_incOverlapData = false);

  static void
  computeL1Norm(Real&                                 a_norm,
                  const Vector<LevelData<EBCellFAB>* >& a_q,
                  const Vector<Real>&                   a_dxVect,
                  const Vector<int>&                    a_refRatio,
                  bool                                  a_incOverlapData = false);


private:

};

#endif

