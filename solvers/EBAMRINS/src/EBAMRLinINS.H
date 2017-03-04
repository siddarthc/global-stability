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

  void defineBaseflowVelocityGradients();

  void defineNewBaseflowVel(const int a_startLevel = -1);

  void setupCoveredBaseAdvVelocity(const Vector<LevelData<EBCellFAB>* >& a_baseVelo,
                                   const Vector<LevelData<EBFluxFAB>* >& a_baseAdvVelo);

  void computeBaseflowVelocityGradients();

//  void advance();

//  void predictor();
//  void corrector();

  virtual void transverseVelocityPredictor(Vector<LevelData<EBCellFAB>* >& a_uGradS,
                                           Vector<LevelData<EBCellFAB>* >& a_scalOld,
                                           bool                            a_reallyVelocity);

  virtual void correctVelocity();

  virtual void computeExtraSourceForPredictor(Vector<LevelData<EBCellFAB>* > & a_source, const Vector<LevelData<EBCellFAB>* >& a_velo, const int& a_dir);

  virtual void computeExtraSourceForCorrector(Vector<LevelData<EBCellFAB>* > & a_source, const Vector<LevelData<EBCellFAB>* >& a_velo);

  // do dot poduct at each point
  // output: a_out(i) = (a_in1(0)*a_in2(0) + a_in1(1)*a_in2(1)+...) at i
  void pointwiseDotProduct(Vector<LevelData<EBCellFAB>* >& a_out,
                           const Vector<LevelData<EBCellFAB>* >& a_in1,
                           const Vector<LevelData<EBCellFAB>* >& a_in2);

  void conclude();

  bool   m_doAdjoint;

  EBIBC* m_baseflowIBC;

  //baseflow data
  Vector<LevelData<EBCellFAB>* > m_baseVelo;
  Vector<LevelData<EBFluxFAB>* > m_baseAdvVelo; // baseflow advective velocity
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredBaseAdvVelLo;
  Vector<LayoutData< Vector< BaseIVFAB<Real> * > >* >  m_coveredBaseAdvVelHi;

  // baseflow velocity gradient
  Vector<Vector<LevelData<EBCellFAB>* > > m_baseVeloGrad;

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
