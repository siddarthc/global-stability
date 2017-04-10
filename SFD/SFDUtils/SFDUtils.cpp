/*
*
*
*
*/

#include "SFDUtils.H"

void SFDUtils::
computeL2Norm(Real&                                 a_norm,
              const Vector<LevelData<EBCellFAB>* >& a_q,
              const Vector<ChomboSFDInterface>&     a_SFDInterface,
              const Vector<Real>&                   a_dxVect,
              const Vector<int>&                    a_refRatio,
              int                                   a_filterIndex,
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
    const LevelData<EBCellFAB>& finestQBar = *((*(a_SFDInterface[finestLevel].getData()))[a_filterIndex]);
    int nvar = finestQ.nComp();
    const DisjointBoxLayout& finestGrids = finestQ.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& Q = finestQ[finestDit()];
      const EBCellFAB& QBar = finestQBar[finestDit()];
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
          ke += (Q(vofit(), ivar) - QBar(vofit(), ivar))*(Q(vofit(), ivar) - QBar(vofit(), ivar));
        }
        
        locNorm += volfrac*cellVol*sqrt(ke);
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
    const LevelData<EBCellFAB>& levelQBar = *((*(a_SFDInterface[ilev].getData()))[a_filterIndex]);
    int nvar = levelQ.nComp();
    const DisjointBoxLayout& levelGrids = levelQ.getBoxes();
    const DisjointBoxLayout& finerGrids = a_q[ilev+1]->getBoxes(); 

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& Q = levelQ[levelDit()];
      const EBCellFAB& QBar = levelQBar[levelDit()];
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
          ke += (Q(vofit(), ivar) - QBar(vofit(), ivar))*(Q(vofit(), ivar) - QBar(vofit(), ivar));
        }

        locNorm += volfrac*cellVol*sqrt(ke);
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

  a_norm /= domainVolume;
}
/**********/

void SFDUtils::
computeLinfNorm(Real&                                 a_norm,
                const Vector<LevelData<EBCellFAB>* >& a_q,
                const Vector<ChomboSFDInterface>&     a_SFDInterface,
                const Vector<Real>&                   a_dxVect,
                const Vector<int>&                    a_refRatio,
                int                                   a_filterIndex,
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
    const LevelData<EBCellFAB>& finestQBar = *((*(a_SFDInterface[finestLevel].getData()))[a_filterIndex]);
    int nvar = finestQ.nComp();
    const DisjointBoxLayout& finestGrids = finestQ.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& Q = finestQ[finestDit()];
      const EBCellFAB& QBar = finestQBar[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = Q.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        for (int ivar = 0; ivar < nComp; ivar++)
        {
          Real val = abs(Q(vofit(), ivar) - QBar(vofit(), ivar));
          if (val > locNorm) locNorm = val;
        }
        
      }
    }
  } // end scoping
 
  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    const LevelData<EBCellFAB>& levelQ = *a_q[ilev];
    const LevelData<EBCellFAB>& levelQBar = *((*(a_SFDInterface[ilev].getData()))[a_filterIndex]);
    int nvar = levelQ.nComp();
    const DisjointBoxLayout& levelGrids = levelQ.getBoxes();
    const DisjointBoxLayout& finerGrids = a_q[ilev+1]->getBoxes(); 

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& Q = levelQ[levelDit()];
      const EBCellFAB& QBar = levelQBar[levelDit()];
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
          Real val = abs(Q(vofit(), ivar) - QBar(vofit(), ivar));
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
