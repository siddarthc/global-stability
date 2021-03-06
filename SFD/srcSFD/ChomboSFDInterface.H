/*
 *
 *
 *
 */


#ifndef _ChomboSFDInterface_H_
#define _ChomboSFDInterface_H_

#include "EBLevelGrid.H"
#include "Vector.H"
#include "EBPWLFineInterp.H"
#include "SFDOp.H"
#include "ChomboSFDOpFunctions.H"
#include "EBCoarseAverage.H"

class ChomboSFDInterface
{

  typedef void (SFDOpFunctions<dataType_>::*base_aXbY_)(dataType_&, const dataType_&, const double&, const dataType_&, const double&);
  typedef void (SFDOpFunctions<dataType_>::*base_copier_)(dataType_&, const dataType_&);

public:

  ChomboSFDInterface();

  ~ChomboSFDInterface();

  void define(int    a_nFilters,
              double a_smallestFilter,
              double a_largestFilter,
              double a_controlCoef,
              int    a_nComp);

  void initialize(const LevelData<EBCellFAB>& a_data,
                  EBLevelGrid&          a_eblg);

  void postInitialize(ChomboSFDInterface* a_finerLevel,
                      EBCoarseAverage& a_ebCoarseAverage);

  void doEBCoarseAverage(ChomboSFDInterface* a_finerLevel,
                         EBCoarseAverage& a_ebCoarseAverage);

  void operator() (LevelData<EBCellFAB>& a_q, double a_dt, EBLevelGrid& a_eblg);

  void regrid(LevelData<EBCellFAB>&                 a_tempFAB,
              ChomboSFDInterface*                   a_coarLevel,
              EBLevelGrid&                          a_eblg,
              EBPWLFineInterp&                      a_interp);

  virtual Vector<string> getDataNames() const;

  virtual Vector<string> getDataNames(int a_filterIndex) const;

  virtual void appendDataToFAB(LevelData<FArrayBox>& a_result, 
                               int                   a_startIndex,
                               const EBLevelGrid&          a_eblg) const;

  virtual void appendDataToFAB(LevelData<FArrayBox>&       a_result, 
                               const LevelData<EBCellFAB>& a_data,
                               int                         a_startIndex, 
                               const EBLevelGrid&                a_eblg) const;

  Vector<LevelData<EBCellFAB>* >* getData();

  bool convergedToSteadyState(const LevelData<EBCellFAB>& a_data, EBLevelGrid& a_eblg);

  int getnComp() const;

// HDF5 stuff
#ifdef CH_USE_HDF5

  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  virtual void writeCheckpointLevel(HDF5Handle& a_handle, int a_filterIndex) const;

  virtual void readCheckpointLevel(HDF5Handle& a_handle, EBLevelGrid& a_eblg);

  virtual void readCheckpointLevel(HDF5Handle& a_handle, int a_filterIndex, EBLevelGrid& a_eblg);

  virtual void writePlotHeader(HDF5HeaderData& a_header, int a_startIndex) const;

  virtual void writePlotHeader(HDF5HeaderData& a_header, int a_comp, int a_startIndex) const;

#endif


protected:

  bool m_isDefined;
  bool m_isInitialized;

  int m_nComp;

  Vector<LevelData<EBCellFAB>* > m_qBar;

  SFDOp<LevelData<EBCellFAB> > m_SFDOp;

};

#endif
