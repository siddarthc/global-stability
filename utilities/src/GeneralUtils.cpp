/*
*
*
*
*/

#include "GeneralUtils.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "EBLevelDataOps.H"
#include "EBAMRDataOps.H"
#include "EBAMRIO.H"
#include <cmath>

/*********************/


void GeneralUtils::
computeLinfNorm(Real&                                 a_norm,
                const Vector<LevelData<EBCellFAB>* >& a_q,
                const Vector<Real>&                   a_dxVect,
                const Vector<int>&                    a_refRatio,
                bool                                  a_incOverlapData)
{
  int finestLevel = a_q.size() - 1;

  a_norm = 0.;
  Real locNorm = 0.; // norm on this processor
  int nComp = a_q[0]->nComp();

  // First, copy the finest level's data
  // scoping begin
  {
    const LevelData<EBCellFAB>& finestQ = *a_q[finestLevel];
    int nvar = finestQ.nComp();
    const DisjointBoxLayout& finestGrids = finestQ.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& Q = finestQ[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = abs(Q(vofit(), ivar));
          if (val > locNorm) locNorm = val;
        }

      }
    }
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {

    const LevelData<EBCellFAB>& levelQ = *a_q[ilev];
    int nvar = levelQ.nComp();
    const DisjointBoxLayout& levelGrids = levelQ.getBoxes();
    const DisjointBoxLayout& finerGrids = a_q[ilev+1]->getBoxes();

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& Q = levelQ[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(thisBox);

      //need to remove the IVs that are already counted at the finerLevel
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

      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = abs(Q(vofit(), ivar) );
          if (val > locNorm) locNorm = val;
        }

      }
    }
  }

  // max among all procs:
  int baseProc = 0;
  Vector<Real> normVec;
  gather(normVec, locNorm, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(normVec.size() == numProc());
    for (int ivec = 0; ivec < normVec.size(); ivec++)
    {
      if (normVec[ivec] > a_norm) a_norm = normVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(a_norm, baseProc);
}

/***********************************************/

/* Compute L0 Norm */


/***********************************************/
void GeneralUtils::
computeL0Norm(Real&                                 a_norm,
              const Vector<LevelData<EBCellFAB>* >& a_q,
              const Vector<Real>&                   a_dxVect,
              const Vector<int>&                    a_refRatio,
              bool                                  a_incOverlapData)
{
  int finestLevel = a_q.size() - 1;

  a_norm = 0.;
  Real locNorm = 0.; // norm on this processor
  Real domainVolume = 0.;
  Real locDomVolume = 0.;
  int nComp = a_q[0]->nComp();

  // First, copy the finest level's data
  // scoping begin
  {
    Real cellVol = a_dxVect[finestLevel];
    for (int idir = 1; idir < SpaceDim; idir++)
    {
      cellVol *= a_dxVect[finestLevel];
    }

    const LevelData<EBCellFAB>& finestQ = *a_q[finestLevel];

    int nvar = finestQ.nComp();
    const DisjointBoxLayout& finestGrids = finestQ.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& Q = finestQ[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)

 {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0;
        for (int ivar = 0; ivar < nComp; ivar++)
        {
           
		if(Q(vofit(), ivar)!=0)
                {

		 ke +=1;// (Q(vofit(), ivar) )*(Q(vofit(), ivar));

		}
        }
        locNorm += volfrac*cellVol*ke;
      }
    }
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    Real cellVol = a_dxVect[ilev];
    for (int idir = 1; idir < SpaceDim; idir++)
    {
      cellVol *= a_dxVect[ilev];
    }

    const LevelData<EBCellFAB>& levelQ = *a_q[ilev];
    int nvar = levelQ.nComp();
    const DisjointBoxLayout& levelGrids = levelQ.getBoxes();
    const DisjointBoxLayout& finerGrids = a_q[ilev+1]->getBoxes();

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& Q = levelQ[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(thisBox);

      //need to remove the IVs that are already counted at the finerLevel
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



 for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0;
        for (int ivar = 0; ivar < nComp; ivar++)
        {
		if(Q(vofit(), ivar)!=0)
		{

	           ke +=1;// (Q(vofit(), ivar));
		}

        }

        locNorm += volfrac*cellVol*ke;
      }
    }
  }

  int baseProc = 0;
  Vector<Real> domVolVec;
  Vector<Real> normVec;
  gather(domVolVec, locDomVolume, baseProc);
  gather(normVec, locNorm, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(domVolVec.size() == numProc());
    for (int ivec = 0; ivec < domVolVec.size(); ivec++)
    {
      domainVolume += domVolVec[ivec];
      a_norm += normVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(domainVolume, baseProc);
  broadcast(a_norm, baseProc);

  a_norm = (a_norm);
  a_norm /= domainVolume;
}


/***********************************************/

/* Compute L1 Norm */


/***********************************************/

void GeneralUtils::
computeL1Norm(Real&                                 a_norm,
              const Vector<LevelData<EBCellFAB>* >& a_q,
              const Vector<Real>&                   a_dxVect,
              const Vector<int>&                    a_refRatio,
              bool                                  a_incOverlapData)
{
  int finestLevel = a_q.size() - 1;

  a_norm = 0.;
  Real locNorm = 0.; // norm on this processor
  Real domainVolume = 0.;
  Real locDomVolume = 0.;
  int nComp = a_q[0]->nComp();

  // First, copy the finest level's data
  // scoping begin
  {
    Real cellVol = a_dxVect[finestLevel];
    for (int idir = 1; idir < SpaceDim; idir++)
    {
      cellVol *= a_dxVect[finestLevel];
    }

    const LevelData<EBCellFAB>& finestQ = *a_q[finestLevel];

    int nvar = finestQ.nComp();
    const DisjointBoxLayout& finestGrids = finestQ.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& Q = finestQ[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0;
        for (int ivar = 0; ivar < nComp; ivar++)
        {
             ke += abs((Q(vofit(), ivar) ));

        }
        locNorm += volfrac*cellVol*ke;
      }
    }
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    Real cellVol = a_dxVect[ilev];
    for (int idir = 1; idir < SpaceDim; idir++)

    {
      cellVol *= a_dxVect[ilev];
    }

    const LevelData<EBCellFAB>& levelQ = *a_q[ilev];
    int nvar = levelQ.nComp();
    const DisjointBoxLayout& levelGrids = levelQ.getBoxes();
    const DisjointBoxLayout& finerGrids = a_q[ilev+1]->getBoxes();

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& Q = levelQ[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(thisBox);

      //need to remove the IVs that are already counted at the finerLevel
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

      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0.;
        for (int ivar = 0; ivar < nComp; ivar++)
        {
           ke += abs((Q(vofit(), ivar)));

        }

        locNorm += volfrac*cellVol*ke;
      }
    }
  }

  int baseProc = 0;
  Vector<Real> domVolVec;
  Vector<Real> normVec;
  gather(domVolVec, locDomVolume, baseProc);
  gather(normVec, locNorm, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(domVolVec.size() == numProc());
    for (int ivec = 0; ivec < domVolVec.size(); ivec++)
    {
      domainVolume += domVolVec[ivec];
      a_norm += normVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(domainVolume, baseProc);
  broadcast(a_norm, baseProc);

  a_norm = (a_norm);
  a_norm /= domainVolume;
}

/***********************************************/

/* Compute L2 Norm */


/***********************************************/


/**********************************************************************/
/**/
void GeneralUtils::
computeL2Norm(Real&                                 a_norm,
              const Vector<LevelData<EBCellFAB>* >& a_q,
              const Vector<Real>&                   a_dxVect,
              const Vector<int>&                    a_refRatio,
              bool                                  a_incOverlapData)
{
  int finestLevel = a_q.size() - 1;

  a_norm = 0.;
  Real locNorm = 0.; // norm on this processor
  Real domainVolume = 0.;
  Real locDomVolume = 0.;
  int nComp = a_q[0]->nComp();

  // First, copy the finest level's data
  // scoping begin
  {
    Real cellVol = a_dxVect[finestLevel];
    for (int idir = 1; idir < SpaceDim; idir++)
    {
      cellVol *= a_dxVect[finestLevel];
    }

    const LevelData<EBCellFAB>& finestQ = *a_q[finestLevel];
    
    int nvar = finestQ.nComp();
    const DisjointBoxLayout& finestGrids = finestQ.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& Q = finestQ[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0;
        for (int ivar = 0; ivar < nComp; ivar++)
        {
             ke += (Q(vofit(), ivar) )*(Q(vofit(), ivar));

        }

        locNorm += volfrac*cellVol*cellVol*ke;
      }
    }
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    Real cellVol = a_dxVect[ilev];
    for (int idir = 1; idir < SpaceDim; idir++)
    {
      cellVol *= a_dxVect[ilev];
    }

    const LevelData<EBCellFAB>& levelQ = *a_q[ilev];
    int nvar = levelQ.nComp();
    const DisjointBoxLayout& levelGrids = levelQ.getBoxes();
    const DisjointBoxLayout& finerGrids = a_q[ilev+1]->getBoxes();

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& Q = levelQ[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(thisBox);

      //need to remove the IVs that are already counted at the finerLevel
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

      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0.;
        for (int ivar = 0; ivar < nComp; ivar++)
        {
    	   ke += (Q(vofit(), ivar) )*(Q(vofit(), ivar));

        }

        locNorm += volfrac*cellVol*cellVol*ke;
      }
    }
  }

  int baseProc = 0;
  Vector<Real> domVolVec;
  Vector<Real> normVec;
  gather(domVolVec, locDomVolume, baseProc);
  gather(normVec, locNorm, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(domVolVec.size() == numProc());
    for (int ivec = 0; ivec < domVolVec.size(); ivec++)
    {
      domainVolume += domVolVec[ivec];
      a_norm += normVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(domainVolume, baseProc);
  broadcast(a_norm, baseProc);

  a_norm = sqrt(a_norm);
  a_norm /= domainVolume;
}

/*********/
void GeneralUtils::
computeKappaDotProduct(Real& a_dotProd,
                       const Vector<LevelData<EBCellFAB>* >& a_in1,
                       const Vector<LevelData<EBCellFAB>* >& a_in2,
                       Vector<int> a_refRatio,
                       Real a_coarsestDx)
{
  a_dotProd = 0.;
  Real locDotProd = 0.; // dot prod on this processor
  Real domainVolume = 0.;
  Real locDomVol = 0.; // domain volume on this processor
  int finestLevel = a_in1.size() - 1;
  int nComp = a_in1[0]->nComp();

  Real levelDx = a_coarsestDx;
  Real coarsestCellVol = 1.;

  for (int ilev = 1; ilev <= finestLevel; ilev++)
  {
    levelDx /= a_refRatio[ilev-1];
  }

  // First, copy the finest level's data
  // scoping begin
  {
    Real cellVol = 1.;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellVol *= levelDx;
      coarsestCellVol *= a_coarsestDx;
    }

    const LevelData<EBCellFAB>& finestData1 = *a_in1[finestLevel];
    const LevelData<EBCellFAB>& finestData2 = *a_in2[finestLevel];
    const DisjointBoxLayout& finestGrids = finestData1.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& data1 = finestData1[finestDit()];
      const EBCellFAB& data2 = finestData2[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = data1.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVol += volfrac*cellVol;

        for (int ivar = 0; ivar < nComp; ivar++)
        {
          locDotProd += volfrac*cellVol*data1(vofit(), ivar)*data2(vofit(), ivar);
        }
      }
    }
  } // end scoping

  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    levelDx *= a_refRatio[ilev];
    if (ilev == 0) CH_assert(levelDx == a_coarsestDx);

    Real cellVol = 1.;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
      cellVol *= levelDx;
    }

    const LevelData<EBCellFAB>& levelData1 = *a_in1[ilev];
    const LevelData<EBCellFAB>& levelData2 = *a_in2[ilev];
    const DisjointBoxLayout& levelGrids = levelData1.getBoxes();
    const DisjointBoxLayout& finerGrids = a_in1[ilev+1]->getBoxes();
    
    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& data1 = levelData1[levelDit()];
      const EBCellFAB& data2 = levelData2[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = data1.getEBISBox();
      IntVectSet ivs(thisBox);

      //need to remove the IVs that are already counted at the finerLevel
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

      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVol += volfrac*cellVol;

        for (int ivar = 0; ivar < nComp; ivar++)
        {
          locDotProd += volfrac*cellVol*data1(vofit(), ivar)*data2(vofit(), ivar);
        }
      }
    }
  }

  int baseProc = 0;
  Vector<Real> domVolVec;
  Vector<Real> dotProdVec;
  gather(domVolVec, locDomVol, baseProc);
  gather(dotProdVec, locDotProd, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(domVolVec.size() == numProc());
    for (int ivec = 0; ivec < domVolVec.size(); ivec++)
    {
      domainVolume += domVolVec[ivec];
      a_dotProd += dotProdVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(domainVolume, baseProc);
  broadcast(a_dotProd, baseProc);
  
  a_dotProd /= domainVolume;
}
