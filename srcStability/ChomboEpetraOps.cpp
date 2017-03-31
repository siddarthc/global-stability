/*
 *
 *
 *
 *
 */

#include "ChomboEpetraOps.H"
#include "SPMDI.H"
#include <math.h>
#include "Epetra_ConfigDefs.h"

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/
/*********/
void ChomboEpetraOps::
rollChomboDataToEpetraVec(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                          Epetra_Vector*                        a_EpetraVec,
                          bool                                  a_incOverlapData,
                          const Vector<int>&                    a_refRatio)
{
  int finestLevel = a_ChomboData.size() - 1;
  int curIndex = -1;
  int success = 0;

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
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ((ivar < nvar) && (success == 0)); ivar++)
        {
          curIndex++;
          success = a_EpetraVec->ReplaceMyValues(1, &(ChomboData(vofit(), ivar)),&curIndex);
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
    const DisjointBoxLayout& finerGrids = a_ChomboData[ilev+1]->getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& ChomboData = levelData[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(thisBox);
      
      // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ((ivar < nvar) && (success == 0)); ivar++)
        {
          curIndex++;
          success = a_EpetraVec->ReplaceMyValues(1, &(ChomboData(vofit(), ivar)),&curIndex);
        }

        if (success != 0)
        {
          string error_msg = "ChomboEpetraOps::rollChomboDataToEpetraVec : copy failed for level = " + SSTR(ilev) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }  
    }
  }

  if (curIndex != a_EpetraVec->MyLength()-1)
  { 
    string error_msg = "ChomboEpetraOps::rollChomboDataToEpetraVec : final index = " + SSTR(curIndex) + " is not consistent with the Epetra_Vec length = " + SSTR(a_EpetraVec->MyLength());

    MayDay::Error(error_msg.c_str());
  } 
}
/*********/
/*********/
void ChomboEpetraOps::
unrollEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                            const Epetra_Vector*            a_EpetraVec,
                            bool                            a_isOverlapDataInc,
                            const Vector<int>&              a_refRatio,
                            double                          a_overlapVal)
{
  int finestLevel = a_ChomboData.size() - 1;
  int curIndex = -1;

  // First, copy the finest level's data
  // scoping begin
  {
    LevelData<EBCellFAB>& finestData = *a_ChomboData[finestLevel];
    int nvar = finestData.nComp(); 
    const DisjointBoxLayout& finestGrids = finestData.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      EBCellFAB& ChomboData = finestData[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          curIndex++;
          ChomboData(vofit(), ivar) = (*a_EpetraVec)[curIndex];
        }
        
      }
    }  
  } // end scoping 

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    LevelData<EBCellFAB>& levelData = *a_ChomboData[ilev];
    int nvar = levelData.nComp();
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    const DisjointBoxLayout& finerGrids = a_ChomboData[ilev+1]->getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      EBCellFAB& ChomboData = levelData[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(thisBox);
      
      // if a_isOverlapDataInc is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_isOverlapDataInc)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, set the value at the overlap IVs to 0
          if (!overlapBox.isEmpty())
          {
            ChomboData.setVal(a_overlapVal, overlapBox, 0, nvar);
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          curIndex++;
          ChomboData(vofit(), ivar) = (*a_EpetraVec)[curIndex];
        }

      }  
    }
  }

  if (curIndex != a_EpetraVec->MyLength()-1)
  { 
    string error_msg = "ChomboEpetraOps::unrollChomboDataToEpetraVec : final index = " + SSTR(curIndex) + " is not consistent with the Epetra_Vec length = " + SSTR(a_EpetraVec->MyLength());

    MayDay::Error(error_msg.c_str());
  }
}
/*********/
/*********/
void ChomboEpetraOps::
scaleEpetraVecByChomboData(Epetra_Vector*                        a_v,
                           const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                           const double&                         a_scaleCD,
                           const int&                            a_startCompV,
                           const int&                            a_startCompCD,
                           const int&                            a_nComp,
                           const int&                            a_totComp,
                           bool                                  a_incOverlapData,
                           const Vector<int>&                    a_refRatio)
{
  CH_assert(a_ChomboData[0]->nComp() >= a_nComp);

  int finestLevel = a_ChomboData.size() - 1;
  int skip = a_totComp - a_nComp;
  int curIndex = a_startCompV - skip - 1;
  int success = 0;

  // First, copy the finest level's data
  // scoping begin
  {
    const LevelData<EBCellFAB>& finestData = *a_ChomboData[finestLevel];
    const DisjointBoxLayout& finestGrids = finestData.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& ChomboData = finestData[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = a_startCompCD; ((ivar < a_nComp) && (success == 0)); ivar++)
        {
          curIndex++;
          double value = (*a_v)[curIndex]; 
          value *= a_scaleCD*ChomboData(vofit(), ivar);
          success = a_v->ReplaceMyValues(1, &value, &curIndex);
        }
        
        if (success != 0) 
        {
          string error_msg = "ChomboEpetraOps::scaleEpetraVecByChomboData : copy failed for level = " + SSTR(finestLevel) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }
    }  
  } // end scoping  

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    const LevelData<EBCellFAB>& levelData = *a_ChomboData[ilev];
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    const DisjointBoxLayout& finerGrids = a_ChomboData[ilev+1]->getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& ChomboData = levelData[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(thisBox);
      
      // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = 0; ((ivar < a_nComp) && (success == 0)); ivar++)
        {
          curIndex++;
          double value = (*a_v)[curIndex];
          value *= a_scaleCD*ChomboData(vofit(), ivar);
          success = a_v->ReplaceMyValues(1, &value, &curIndex);
        }

        if (success != 0)
        {
          string error_msg = "ChomboEpetraOps::scaleEpetraVecByChomboData : copy failed for level = " + SSTR(ilev) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }  
    }
  } 

  int checker = curIndex + a_totComp - a_startCompV - a_nComp;

  if (checker != a_v->MyLength()-1)
  { 
    string error_msg = "ChomboEpetraOps::scaleEpetraVecByChomboData : final index = " + SSTR(checker) + " is not consistent with the Epetra_Vec length = " + SSTR(a_v->MyLength());

    MayDay::Error(error_msg.c_str());
  } 
}
/*********/
/*********/
void ChomboEpetraOps::
divideEpetraVecByChomboData(Epetra_Vector*                        a_v,
                            const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                            const double&                         a_scaleCD,
                            const int&                            a_startCompV,
                            const int&                            a_startCompCD,
                            const int&                            a_nComp,
                            const int&                            a_totComp,
                            bool                                  a_incOverlapData,
                            const Vector<int>&                    a_refRatio)
{
  CH_assert(a_ChomboData[0]->nComp() >= a_nComp);

  int finestLevel = a_ChomboData.size() - 1;
  int skip = a_totComp - a_nComp;
  int curIndex = a_startCompV - skip - 1;
  int success = 0;

  // First, copy the finest level's data
  // scoping begin
  {
    const LevelData<EBCellFAB>& finestData = *a_ChomboData[finestLevel];
    const DisjointBoxLayout& finestGrids = finestData.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& ChomboData = finestData[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = a_startCompCD; ((ivar < a_nComp) && (success == 0)); ivar++)
        {
          curIndex++;
          double value = (*a_v)[curIndex]; 
          double denom = (a_scaleCD*ChomboData(vofit(), ivar));
          value *= 1./denom;
          success = a_v->ReplaceMyValues(1, &value, &curIndex);
        }
        
        if (success != 0) 
        {
          string error_msg = "ChomboEpetraOps::divideEpetraVecByChomboData : copy failed for level = " + SSTR(finestLevel) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }
    }  
  } // end scoping  

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    const LevelData<EBCellFAB>& levelData = *a_ChomboData[ilev];
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    const DisjointBoxLayout& finerGrids = a_ChomboData[ilev+1]->getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& ChomboData = levelData[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(thisBox);
      
      // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = 0; ((ivar < a_nComp) && (success == 0)); ivar++)
        {
          curIndex++;
          double value = (*a_v)[curIndex];
          double denom = (a_scaleCD*ChomboData(vofit(), ivar));
          value *= 1./denom;
          success = a_v->ReplaceMyValues(1, &value, &curIndex);
        }

        if (success != 0)
        {
          string error_msg = "ChomboEpetraOps::divideEpetraVecByChomboData : copy failed for level = " + SSTR(ilev) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }  
    }
  } 

  int checker = curIndex + a_totComp - a_startCompV - a_nComp;

  if (checker != a_v->MyLength()-1)
  { 
    string error_msg = "ChomboEpetraOps::divideEpetraVecByChomboData : final index = " + SSTR(checker) + " is not consistent with the Epetra_Vec length = " + SSTR(a_v->MyLength());

    MayDay::Error(error_msg.c_str());
  } 
}
/*********/
/*********/
void ChomboEpetraOps::
addChomboDataToEpetraVec(Epetra_Vector*                        a_v,
                         const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                         const double&                         a_scaleV,
                         const double&                         a_scaleCD,
                         const int&                            a_startCompV,
                         const int&                            a_startCompCD,
                         const int&                            a_nComp,
                         const int&                            a_totComp,
                         bool                                  a_incOverlapData,
                         const Vector<int>&                    a_refRatio)
{
  CH_assert(a_ChomboData[0]->nComp() >= a_nComp);

  int finestLevel = a_ChomboData.size() - 1;
  int skip = a_totComp - a_nComp;
  int curIndex = a_startCompV - skip - 1;
  int success = 0;

  // First, copy the finest level's data
  // scoping begin
  {
    const LevelData<EBCellFAB>& finestData = *a_ChomboData[finestLevel];
    const DisjointBoxLayout& finestGrids = finestData.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& ChomboData = finestData[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = a_startCompCD; ((ivar < a_nComp) && (success == 0)); ivar++)
        {
          curIndex++;
          double value = a_scaleV*(*a_v)[curIndex]; 
          value += a_scaleCD*ChomboData(vofit(), ivar);
          success = a_v->ReplaceMyValues(1, &value, &curIndex);
        }
        
        if (success != 0) 
        {
          string error_msg = "ChomboEpetraOps::addChomboDataToEpetraVec : copy failed for level = " + SSTR(finestLevel) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }
    }  
  } // end scoping  

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    const LevelData<EBCellFAB>& levelData = *a_ChomboData[ilev];
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    const DisjointBoxLayout& finerGrids = a_ChomboData[ilev+1]->getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& ChomboData = levelData[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(thisBox);
      
      // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = 0; ((ivar < a_nComp) && (success == 0)); ivar++)
        {
          curIndex++;
          double value = a_scaleV*(*a_v)[curIndex];
          value += a_scaleCD*ChomboData(vofit(), ivar);
          success = a_v->ReplaceMyValues(1, &value, &curIndex);
        }

        if (success != 0)
        {
          string error_msg = "ChomboEpetraOps::addChomboDataToEpetraVec : copy failed for level = " + SSTR(ilev) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }  
    }
  } 

  int checker = curIndex + a_totComp - a_startCompV - a_nComp;

  if (checker != a_v->MyLength()-1)
  { 
    string error_msg = "ChomboEpetraOps::addChomboDataToEpetraVec : final index = " + SSTR(checker) + " is not consistent with the Epetra_Vec length = " + SSTR(a_v->MyLength());

    MayDay::Error(error_msg.c_str());
  } 
}
/*********/
/*********/
void ChomboEpetraOps::
addEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* > & a_ChomboData,
                         const Epetra_Vector*             a_v,
                         const double&                    a_scaleCD,
                         const double&                    a_scaleV,
                         const int&                       a_startCompCD,
                         const int&                       a_startCompV,
                         const int&                       a_nComp,
                         const int&                       a_totComp,
                         bool                             a_incOverlapData,
                         const Vector<int>&               a_refRatio)
{
  CH_assert(a_ChomboData[0]->nComp() >= a_nComp);

  int finestLevel = a_ChomboData.size() - 1;
  int skip = a_totComp - a_nComp;
  int curIndex = a_startCompV - skip - 1;

  // First, copy the finest level's data
  // scoping begin
  {
    LevelData<EBCellFAB>& finestData = *a_ChomboData[finestLevel];
    const DisjointBoxLayout& finestGrids = finestData.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      EBCellFAB& ChomboData = finestData[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = a_startCompCD; ivar < a_nComp; ivar++)
        {
          curIndex++;
          double value = a_scaleV*(*a_v)[curIndex];
          value += a_scaleCD*ChomboData(vofit(), ivar);
          ChomboData(vofit(), ivar) = value;
        }
        
      }
    }  
  } // end scoping 
 
  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    LevelData<EBCellFAB>& levelData = *a_ChomboData[ilev];
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    const DisjointBoxLayout& finerGrids = a_ChomboData[ilev+1]->getBoxes();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      EBCellFAB& ChomboData = levelData[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = ChomboData.getEBISBox();
      IntVectSet ivs(thisBox);
      
      // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        curIndex += skip;
        for (int ivar = 0; ivar < a_nComp; ivar++)
        {
          curIndex++;
          double value = a_scaleV*(*a_v)[curIndex];
          value += a_scaleCD*ChomboData(vofit(), ivar);
          ChomboData(vofit(), ivar) = value;
        }

      }  
    }
  } 

  int checker = curIndex + a_totComp - a_startCompV - a_nComp;
  if (checker != a_v->MyLength()-1)
  { 
    string error_msg = "ChomboEpetraOps::addEpetraVecToChomboData : final index = " + SSTR(checker) + " is not consistent with the Epetra_Vec length = " + SSTR(a_v->MyLength());

    MayDay::Error(error_msg.c_str());
  } 
}
/*********/
/*********/
int ChomboEpetraOps::
getnElementsOnThisProc(const Vector<DisjointBoxLayout>& a_ChomboDBL,
                       const Vector<EBLevelGrid>&       a_ChomboEBLG,
                       const int&                       a_nComp,
                       bool                             a_incOverlapData,
                       const Vector<int>&               a_refRatio)
{
  int nLocElem = 0;
  int finestLevel = a_ChomboDBL.size() - 1;

  // First, copy the finest level's data
  // scoping begin
  {
    int nvar = a_nComp; 
    const DisjointBoxLayout& finestGrids = a_ChomboDBL[finestLevel];
    const EBISLayout&        finestEBISL = a_ChomboEBLG[finestLevel].getEBISL();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = finestEBISL[finestDit()];
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          nLocElem++;
        }
        
      }
    }  
  } // end scoping 

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    int nvar = a_nComp;
    const DisjointBoxLayout& levelGrids = a_ChomboDBL[ilev];
    const DisjointBoxLayout& finerGrids = a_ChomboDBL[ilev+1];
    const EBISLayout&        levelEBISL = a_ChomboEBLG[ilev].getEBISL();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = levelEBISL[levelDit()];
      IntVectSet ivs(thisBox);
      
     // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          nLocElem++;
        }

      }  
    }
  } 

  return nLocElem;
}
/*********/
/*********/
Epetra_Map ChomboEpetraOps::
getEpetraMap(const Vector<DisjointBoxLayout>& a_ChomboDBL,
             const Vector<EBLevelGrid>&       a_ChomboEBLG,
             const int&                       a_nComp,
             const Epetra_Comm*               a_commPtr,
             bool                             a_incOverlapData,
             const Vector<int>&               a_refRatio)
{
  int nLocElem = getnElementsOnThisProc(a_ChomboDBL, a_ChomboEBLG, a_nComp, a_incOverlapData, a_refRatio);
  Epetra_Map retval(-1, nLocElem, 0, *a_commPtr);
  return retval;
}
/*********/
/*********/
void ChomboEpetraOps::
getVolWeights(Epetra_Vector*                   a_weights,
              double&                          a_domainVolume,
              const Vector<DisjointBoxLayout>& a_ChomboDBL,
              const Vector<EBLevelGrid>&       a_ChomboEBLG,
              const int&                       a_nComp,
              const double&                    a_coarsestDx,
              const Vector<int>&               a_refRatio,
              bool                             a_incOverlapData)
{
  int finestLevel = a_ChomboDBL.size() - 1;
  int curIndex = -1;
  int success = 0;

  Real levelDx = a_coarsestDx;
  Real coarsestCellVol = 1.;

  a_domainVolume = 0.;  
  Real locDomVol = 0.; // domain volume on this processor

  // First, copy the finest level's data
  // scoping begin
  {
    for (int ilev = 1; ilev <= finestLevel; ilev++)
    {
      levelDx /= a_refRatio[ilev-1];
    }

    Real cellvol = 1.;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellvol *= levelDx;
      coarsestCellVol *= a_coarsestDx;
    }

    int nvar = a_nComp; 
    const DisjointBoxLayout& finestGrids = a_ChomboDBL[finestLevel];
    const EBISLayout&        finestEBISL = a_ChomboEBLG[finestLevel].getEBISL();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = finestEBISL[finestDit()];
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        Real weight = sqrt(volfrac*cellvol);
        locDomVol += volfrac*cellvol;

        for (int ivar = 0; ((ivar < nvar) && (success == 0)); ivar++)
        {
          curIndex++;
          success = a_weights->ReplaceMyValues(1, &weight ,&curIndex);
        }
        
        if (success != 0) 
        {
          string error_msg = "ChomboEpetraOps::getVolWeights : copy failed for level = " + SSTR(finestLevel) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }
    }  
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    levelDx *= a_refRatio[ilev];
    if (ilev == 0) CH_assert(levelDx == a_coarsestDx);

    Real cellvol = 1.;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellvol *= levelDx;
    }

    int nvar = a_nComp;
    const DisjointBoxLayout& levelGrids = a_ChomboDBL[ilev];
    const DisjointBoxLayout& finerGrids = a_ChomboDBL[ilev+1];
    const EBISLayout&        levelEBISL = a_ChomboEBLG[ilev].getEBISL();
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = levelEBISL[levelDit()];
      IntVectSet ivs(thisBox);
      
      // if a_incOverlapData is false, need to remove the IVs that are already counted at the finerLevel
      if (!a_incOverlapData)
      {
        for (finerLit.begin(); finerLit.ok(); ++finerLit)
        {
          Box overlapBox = finerGrids[finerLit];
          overlapBox.coarsen(a_refRatio[ilev]);
          overlapBox &= thisBox;
          // if overlap, remove the overlap IVs
          if (!overlapBox.isEmpty())
          {
            IntVectSet ivsExclude(overlapBox);
            ivs -= ivsExclude;
          }
        }
      }
    
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        Real weight = sqrt(volfrac*cellvol);
        locDomVol += volfrac*cellvol;

        for (int ivar = 0; ((ivar < nvar) && (success == 0)); ivar++)
        {
          curIndex++;
          success = a_weights->ReplaceMyValues(1, &weight,&curIndex);
        }

        if (success != 0)
        {
          string error_msg = "ChomboEpetraOps::rollChomboDataToEpetraVec : copy failed for level = " + SSTR(ilev) + " at index = " + SSTR(curIndex);
          MayDay::Error(error_msg.c_str());
        }
      }  
    }
  } 

  int baseProc = 0;
  Vector<Real> domVolVec;
  gather(domVolVec, locDomVol, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(domVolVec.size() == numProc());
    for (int ivec = 0; ivec < domVolVec.size(); ivec++)
    {
      a_domainVolume += domVolVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(a_domainVolume, baseProc); 

  double scalingFactor = sqrt(1./coarsestCellVol);
  a_weights->Scale(scalingFactor);
}
/*********/
/*********/
int ChomboEpetraOps::
computeL2Norm(double&              a_norm,
              const Epetra_Vector& a_v)
{
  return a_v.Norm2(&a_norm);
}
/*********/
/*********/
int ChomboEpetraOps::
computeWeightedL2Norm(double&              a_norm,
                      const Epetra_Vector& a_v,
                      const Epetra_Vector& a_weights)
{
  Epetra_Vector tmpVec(a_v);
  tmpVec.PutScalar(0.);
  tmpVec.Multiply(1., a_v, a_weights, 0.); // make tmpVec = a_v*a_weights
  return tmpVec.Norm2(&a_norm);
}
/*********/
/*********/
int ChomboEpetraOps::
computeDotProduct(double&              a_result,
                  const Epetra_Vector& a_v1,
                  const Epetra_Vector& a_v2)
{
  return a_v1.Dot(a_v2, &a_result);
}
/*********/
/*********/
int ChomboEpetraOps::
computeWeightedDotProduct(double& a_result,
                          const Epetra_Vector& a_v1,
                          const Epetra_Vector& a_v2,
                          const Epetra_Vector& a_weights)
{
  Epetra_Vector tmpVec1(a_v1);
  Epetra_Vector tmpVec2(a_v2);
  tmpVec1.PutScalar(0.);
  tmpVec2.PutScalar(0.);
  tmpVec1.Multiply(1., a_v1, a_weights, 0.); // make tmpVec = a_v*a_weights
  tmpVec2.Multiply(1., a_v2, a_weights, 0.); // make tmpVec = a_v*a_weights
  return tmpVec1.Dot(tmpVec2, &a_result);
}
