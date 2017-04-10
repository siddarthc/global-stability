/*
*
*
*
*/

#ifndef _SFD_UTILS_H_
#define _SFD_UTILS_H_

#include "ChomboSFDInterface.H"

class SFDUtils
{
public:

  static void
  computeL2Norm(Real&                                 a_norm,
                const Vector<LevelData<EBCellFAB>* >& a_q,
                const Vector<ChomboSFDInterface>&     a_SFDInterface,
                const Vector<Real>&                   a_dxVect,
                const Vector<int>&                    a_refRatio,
                int                                   a_filterIndex,
                bool                                  a_incOverlapData = false);

  static void
  computeLinfNorm(Real&                                 a_norm,
                  const Vector<LevelData<EBCellFAB>* >& a_q,
                  const Vector<ChomboSFDInterface>&     a_SFDInterface,
                  const Vector<Real>&                   a_dxVect,
                  const Vector<int>&                    a_refRatio,
                  int                                   a_filterIndex,
                  bool                                  a_incOverlapData = false);


private: 

};

#endif
