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


  static void 
  plotSFDIntegratorError(std::string                       a_fileName,
                         const Vector<ChomboSFDInterface>& a_SFDInterface,
                         const Vector<DisjointBoxLayout>&  a_grids,
                         const Vector<EBISLayout>&         a_ebisl,
                         const Vector<ProblemDomain>&      a_domain,
                         const Vector<Real>&               a_dx,
                         const Real&                       a_dt,
                         const Real&                       a_time,
                         const Vector<int>&               a_refRatio,
                         int                               a_filterIndex);

  // this returns m_qdiffNew.m_qdiffOld if PI control is on
  static void
  computeQDiffDotProd(Real&                                 a_dotProd,
      	              const Vector<ChomboSFDInterface>&     a_SFDInterface,
   	              const Vector<Real>&                   a_dxVect,
                      const Vector<int>&                    a_refRatio,
    	              int                                   a_filterIndex,
         	      bool                                  a_incOverlapData = false);

  // this returns L2 norm of integrated error if PI control is on
  static void
  computeIntegratedErrorL2Norm(Real&                                 a_norm,
                               const Vector<ChomboSFDInterface>&     a_SFDInterface,
	                       const Vector<Real>&                   a_dxVect,
	                       const Vector<int>&                    a_refRatio,
	                       int                                   a_filterIndex,
	                       bool                                  a_incOverlapData = false); 

private: 

};

#endif
