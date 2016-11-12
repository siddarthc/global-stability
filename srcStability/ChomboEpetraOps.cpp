/*
 *
 *
 *
 *
 */

#include "ChomboEpetraOps.H"
#include <math.h>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/
/*********/
void ChomboEpetraOps::
rollChomboDataToEpetraVec(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                          Epetra_Vector*                        a_EpetraVec,
                          bool                                  a_incOverlapData,
                          const Vector<Real>&                   a_refRatio)
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
        for (int ivar = 0; (ivar < nvar && success == 0); ivar++)
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
        for (int ivar = 0; (ivar < nvar && success == 0); ivar++)
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
}
/*********/
/*********/
void ChomboEpetraOps::
unrollEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                            const Epetra_Vector*            a_EpetraVec,
                            bool                            a_isOverlapDataInc,
                            const Vector<Real>&             a_refRatio,
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
int ChomboEpetraOps::
getnElementsOnThisProc(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                       bool                                  a_incOverlapData,
                       const Vector<Real>&                   a_refRatio)
{
  int nLocElem = 0;
  int finestLevel = a_ChomboData.size() - 1;

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
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          nLocElem++;
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
getEpetraMap(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
             const Epetra_Comm*                    a_commPtr,
             bool                                  a_incOverlapData,
             const Vector<Real>&                   a_refRatio)
{
  int nLocElem = getnElementsOnThisProc(a_ChomboData, a_incOverlapData);
  Epetra_Map retval(-1, nLocElem, 0, *a_commPtr);
  return retval;
}
/*********/
/*********/
void ChomboEpetraOps::
getVolWeights(Epetra_Vector*                       a_weights,
              double&                              a_domainVolume,
              const Vector<LevelData<EBCellFAB>* > a_ChomboData,
              const double&                        a_coarsestDx,
              const Vector<Real>&                  a_refRatio,
              bool                                 a_incOverlapData)
{
  int finestLevel = a_ChomboData.size() - 1;
  int curIndex = -1;
  int success = 0;

  Real levelDx = a_coarsestDx;

  a_domainVolume = 0.;  

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
    }

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
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        Real weight = sqrt(1./(volfrac*cellvol));
        a_domainVolume += volfrac*cellvol;

        for (int ivar = 0; (ivar < nvar && success == 0); ivar++)
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
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        Real weight = sqrt(1./(volfrac*cellvol));
        a_domainVolume += volfrac*cellvol;

        for (int ivar = 0; (ivar < nvar && success == 0); ivar++)
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

  double scalingFactor = sqrt(a_domainVolume/a_weights->GlobalLength());
  a_weights->Scale(scalingFactor);
}
/*********/
/*********/
int ChomboEpetraOps::
computeL2Norm(double&              a_norm,
              const Epetra_Vector& a_mv)
{
  return a_mv.Norm2(&a_norm);
}
