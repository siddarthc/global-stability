/*
*
*
*
*/

#ifndef _WAVEMAKER_UTILS_H_
#define _WAVEMAKER_UTILS_H_

#include "ProblemDomain.H"
#include "EBIndexSpace.H"
#include "Vector.H"
#include "REAL.H"

// finds the wavemaker given direct and adjoint modes
class WavemakerUtils
{
public:

  static void
  plotWavemaker(std::string               a_fileName,
                std::string               a_directModeFile,
                std::string               a_adjointModeFile,
                Vector<int>               a_refRatio,
                const ProblemDomain&      a_coarsestDomain,
                Real                      a_coarsestDx,
                const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());


  static void 
  computeStructuralSensitivity(Vector<LevelData<EBCellFAB>* >&       a_out,
                               const Vector<LevelData<EBCellFAB>* >& a_directModeVelo,
                               const Vector<LevelData<EBCellFAB>* >& a_adjointModeVelo,
                               Vector<int> a_refRatio,
                               Real a_coarsestDx);

  static void 
  computeKappaDotProduct(Real& a_out,
                         const Vector<LevelData<EBCellFAB>* >& a_in1,
                         const Vector<LevelData<EBCellFAB>* >& a_in2,
                         Vector<int> a_refRatio,
                         Real a_coarsestDx);

private:

/*
  // note that this defines the data
  void getVelocityFromFile(Vector<LevelData<EBCellFAB>* > a_velocity,
                           std::string                    a_fileName,
                           Vector<int>                    a_refRatio,
                           const ProblemDomain&           a_coarsestDomain,
                           const EBIndexSpace* const      a_ebisPtr = Chombo_EBIS::instance());
*/
};

#endif
