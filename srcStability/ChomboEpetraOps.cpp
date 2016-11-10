/*
 *
 *
 *
 *
 */

#include "ChomboEpetraOps.H"

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/
void ChomboEpetraOps::
rollChomboDataToEpetraVec(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                          Epetra_Vector*                        a_EpetraVec,
                          bool                                  a_incOverlapData,
                          const Vector<Real>&                   a_refRatio)
{
  int finestLevel = a_ChomboData.size() - 1;
  int curIndex = 0;

  // First, copy the finest level's data
  // scoping begin
  {
    const LevelData<EBCellFAB>& finestData = *a_ChomboData[finestLevel];
    int nvar = finestData.nComp(); 
    const DisjointBoxLayout& finestGrids = finestData.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& ChomboData = finestData[finestDit()];
      Box box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        int success = 0;
        for (int ivar = 0; (ivar < nvar && success == 0); ivar++)
        {
          success = a_EpetraVec->ReplaceMyValues(1, &(ChomboData(vofit(), ivar)),&curIndex);
          curIndex++;
        }
        
        if (success != 0) 
        {
          string error_msg = "ChomboEpetraOps::rollChomboDataToEpetraVec : copy failed for level = " + SSTR(finestLevel) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }
    }  
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    const LevelData<EBCellFAB>& levelData = *a_ChomboData[ilev];
    int nvar = levelData.nComp();
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
    
    }

  } 
}
/*********/
/*********/
void ChomboEpetraOps::
unrollEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                            const Epetra_Vector*            a_EpetraVec,
                            bool                            a_isOverlapDataInc,
                            const Vector<Real>&             a_refRatio)
{

}
/*********/
/*********/
int ChomboEpetraOps::
getnElementsOnThisProc(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                       bool                                  a_incOverlapData,
                       const Vector<Real>&                   a_refRatio)
{
  return 0;
}
/*********/
/*********/
Epetra_Map ChomboEpetraOps::
getEpetraMap(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
             const Epetra_Comm*                    a_commPtr,
             bool                                  a_incOverlapData,
             const Vector<Real>&                   a_refRatio)
{
  int nLocElem = getnElementsOnThisProc(a_ChomboData, a_incOverlapData);
  Epetra_Map retval(-1, nLocElem, 0, *a_commPtr);
  return retval;
}
