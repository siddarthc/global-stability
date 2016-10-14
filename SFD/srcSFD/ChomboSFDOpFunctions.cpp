/*
 *
 *
 *
 *
 */

#include "ChomboSFDOpFunctions.H"
#include "EBLevelDataOps.H"

ChomboSFDOpFunctions::ChomboSFDOpFunctions()
{
  m_eblgPtr = NULL;
}
/*******/

ChomboSFDOpFunctions::~ChomboSFDOpFunctions()
{
  if (m_eblgPtr != NULL) m_eblgPtr = NULL;
}
/*******/

void ChomboSFDOpFunctions::
aXbY(dataType_&            a_result,
     const dataType_&      a_X,
     const double&         a_aCoef,
     const dataType_&      a_Y,
     const double&         a_bCoef)
{
//  Real aCoef = a_aCoef;
// Real bCoef = a_bCoef;
  EBLevelDataOps::axby(a_result, a_X, a_Y, a_aCoef, a_bCoef);
}
/*******/
void ChomboSFDOpFunctions::
copierFunc(dataType_& a_result,
           const dataType_& a_data)
{
  CH_assert(m_eblgPtr != NULL)

  EBCellFactory fact(m_eblgPtr->getEBISL());
  IntVect ivGhost = m_eblgPtr->getGhost()*IntVect::Unit;
  a_result.define(m_eblgPtr->getDBL(), a_data.nComp(), ivGhost, fact); 

  a_data.copyTo(a_data.interval(), a_result, a_result.interval());
}
/*********/
void ChomboSFDOpFunctions::
setEBLG(EBLevelGrid* a_eblg)
{
  m_eblgPtr = a_eblg;
}
