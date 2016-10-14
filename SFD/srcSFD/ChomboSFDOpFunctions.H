/*
 *
 *
 *
 *
 */

#ifndef _ChomboSFDOpFunctions_H_
#define _ChomboSFDOpFunctions_H_

#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "SFDOpFunctions.H"

typedef LevelData<EBCellFAB> dataType_;

class ChomboSFDOpFunctions : public SFDOpFunctions<dataType_>
{

public:

  ChomboSFDOpFunctions();

  ~ChomboSFDOpFunctions();

  virtual void aXbY(dataType_&            a_result,
                    const dataType_&      a_X,
                    const double&         a_aCoef,
                    const dataType_&      a_Y,
                    const double&         a_bCoef);

  virtual void copierFunc(dataType_&       a_clone,
                          const dataType_& a_data); 


  void setEBLG(EBLevelGrid* a_eblg);

protected:

  EBLevelGrid* m_eblgPtr;

};

#endif
