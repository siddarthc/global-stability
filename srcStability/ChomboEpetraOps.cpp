/*
 *
 *
 *
 *
 */

#include "ChomboEpetraOps.H"

/*********/
void ChomboEpetraOps::
rollChomboDataToEpetraVec(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                          Epetra_Vector*                        a_EpetraVec,
                          bool                                  a_incOverlapData)
{

}
/*********/
/*********/
void ChomboEpetraOps::
unrollEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                            const Epetra_Vector*            a_EpetraVec,
                            bool                            a_isOverlapDataInc)
{

}
/*********/
/*********/
int ChomboEpetraOps::
getnElementsOnThisProc(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                       bool a_incOverlapData)
{
  return 0;
}
/*********/
/*********/
Epetra_Map ChomboEpetraOps::
getEpetraMap(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
             const Epetra_Comm*                    a_commPtr,
             bool a_incOverlapData)
{
  int nLocElem = getnElementsOnThisProc(a_ChomboData, a_incOverlapData);
  Epetra_Map retval(-1, nLocElem, 0, *a_commPtr);
  return retval;
}
