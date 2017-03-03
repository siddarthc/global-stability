/*
 *
 *
 *
 *
 */

#ifndef _EBAMRLININS_H_
#define _EBAMRLININS_H_

#include "EBAMRNoSubcycle.H"

class EBAMRLinINS : public EBAMRNoSubcycle
{

public:

  EBAMRLinINS()
  {
    MayDay::Error("EBAMRLinINS::Invalid Constructor");
  }

  EBAMRLinINS(const AMRParameters&      a_params,
              const EBIBCFactory&       a_baseflowIBC,
              const EBIBCFactory&       a_solverIBC,
              const ProblemDomain&      a_coarsestDomain,
              Real                      a_viscosity,
              bool                      a_doAdjoint,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());


  virtual ~EBAMRLinINS();

  void setupForStabilityRun(const Epetra_Vector&             a_x,
                            const Vector<DisjointBoxLayout>& a_baseflowDBL,
                            const Vector<EBLevelGrid>&       a_baseflowEBLG,
                            const std::string&               a_baseflowFile,
                            double                           a_pertScale,
                            bool                             a_incOverlapData,
                            bool                             a_setupForPlottingData = false);

protected:

  void allocateBaseflowDataHolders();

  void defineBaseflowIrregularData();

  void defineNewBaseflowVel(const int a_startLevel = -1);

  void setupCoveredBaseAdvVelocity(const Vector<LevelData<EBCellFAB>* >& a_baseVelo,
                                   const Vector<LevelData<EBFluxFAB>* >& a_baseAdvVelo);

//  void advance();

//  void predictor();
//  void corrector();

  void conclude();

  bool   m_doAdjoint;

  EBIBC* m_baseflowIBC;

  //baseflow data
  Vector<LevelData<EBCellFAB>* > m_baseVelo;
  Vector<LevelData<EBFluxFAB>* > m_baseAdvVelo; // baseflow advective velocity
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredBaseAdvVelLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredBaseAdvVelHi;

private:

  EBAMRLinINS(const EBAMRLinINS& a_input)
  {
    MayDay::Error("Invalid Constructor");
  }

  void operator=(const EBAMRLinINS& a_input)
  {
    MayDay::Error("Invalid Constructor");
  }
};

#endif
