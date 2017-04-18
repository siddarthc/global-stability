/*
*
*
*
*/

#include "SFDUtils.H"
#include "EBAMRIO.H"
#include "EBCellFactory.H"
#include "EBLevelDataOps.H"

/*********/
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
/*********/
void SFDUtils::
plotSFDIntegratorError(std::string                       a_fileName,
                       const Vector<ChomboSFDInterface>& a_SFDInterface,
                       const Vector<DisjointBoxLayout>&  a_grids,
                       const Vector<EBISLayout>&         a_ebisl,
                       const Vector<ProblemDomain>&      a_domain,
                       const Vector<Real>&               a_dx,
                       const Real&                       a_dt,
                       const Real&                       a_time,
                       const Vector<int>&               a_refRatio,
                       int                               a_filterIndex)
{
#ifdef CH_USE_HDF5

  Vector<LevelData<EBCellFAB>* > tempData;
  tempData.resize(a_SFDInterface.size());
  for (int ilev = 0; ilev < a_SFDInterface.size(); ilev++)
  {
    if (a_SFDInterface[ilev].getIntegratedError(a_filterIndex) == NULL) MayDay::Error("SFDInterface::getIntegratedError returned NULL");
    EBCellFactory ebcellfact(a_ebisl[ilev]);
    tempData[ilev] = new LevelData<EBCellFAB>(a_grids[ilev], SpaceDim,  3*IntVect::Unit, ebcellfact); 
    EBLevelDataOps::assign(*tempData[ilev], *a_SFDInterface[ilev].getIntegratedError(a_filterIndex));
  }

  int curNumLevels = a_SFDInterface.size();
  Vector<string> names(SpaceDim);
  for (int idir = 0; idir < SpaceDim; idir++)
  {
    char velochar[100];
    sprintf(velochar, "integratedError_velo%d", idir);
    names[idir] = string(velochar);
   }

  Vector<Real> coveredValues;
  bool replaceCovered = false;


  writeEBHDF5(a_fileName,
              a_grids,
              tempData,
              names,
              a_domain[0].domainBox(),
              a_dx[0],
              a_dt,
              a_time,
              a_refRatio,
              curNumLevels,
              replaceCovered,
              coveredValues);


//  writeEBAMRname(&tempData, a_fileName);

  for (int ilev = 0; ilev < a_SFDInterface.size(); ilev++)
  {
    if (tempData[ilev] != NULL) delete tempData[ilev];
  }
#else
  MayDay::Error("SFDUtils::plotSFDIntegratorError required HDF5");
#endif
}
/*********/
void SFDUtils::
computeQDiffDotProd(Real&                                 a_dotProd,
                    const Vector<ChomboSFDInterface>&     a_SFDInterface,
                    const Vector<Real>&                   a_dxVect,
                    const Vector<int>&                    a_refRatio,
                    int                                   a_filterIndex,
                    bool                                  a_incOverlapData)
{

  int finestLevel = a_SFDInterface.size() - 1;

  a_dotProd = 0.;
  Real locDotProd = 0.; // norm on this processor

  // First, copy the finest level's data
  // scoping begin
  {
    const LevelData<EBCellFAB>& finestQdiffNew = *(a_SFDInterface[finestLevel].getQdiffNew(a_filterIndex));
    const LevelData<EBCellFAB>& finestQdiffOld = *(a_SFDInterface[finestLevel].getQdiffOld(a_filterIndex));

//    CH_assert((finestQdiffNew != NULL) && (finestQdiffOld != NULL));

    int nvar = finestQdiffNew.nComp();
    const DisjointBoxLayout& finestGrids = finestQdiffNew.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& diffNew = finestQdiffNew[finestDit()];
      const EBCellFAB& diffOld = finestQdiffOld[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = diffNew.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
  
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          locDotProd += (diffNew(vofit(), ivar))*(diffOld(vofit(), ivar));
        }
        
      }
    }
  } // end scoping
 
  for (int ilev = finestLevel-1; ilev >= 0; ilev--)
  {
    const LevelData<EBCellFAB>& levelQdiffNew = *(a_SFDInterface[ilev].getQdiffNew(a_filterIndex));
    const LevelData<EBCellFAB>& levelQdiffOld = *(a_SFDInterface[ilev].getQdiffOld(a_filterIndex));

//    CH_assert((levelQdiffNew != NULL) && (levelQdiffOld != NULL))

    int nvar = levelQdiffNew.nComp();
    const DisjointBoxLayout& levelGrids = levelQdiffNew.getBoxes();
    const DisjointBoxLayout& finerGrids = a_SFDInterface[ilev+1].getQdiffNew(a_filterIndex)->getBoxes(); 

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& diffNew = levelQdiffNew[levelDit()];
      const EBCellFAB& diffOld = levelQdiffOld[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = diffNew.getEBISBox();
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

        for (int ivar = 0; ivar < nvar; ivar++)
        {
          locDotProd += (diffNew(vofit(), ivar))*(diffOld(vofit(), ivar));
        }
      }
    }
  }

  int baseProc = 0;
  Vector<Real> dotVec;
  gather(dotVec, locDotProd, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(dotVec.size() == numProc());
    for (int ivec = 0; ivec < dotVec.size(); ivec++)
    {
      a_dotProd += dotVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(a_dotProd, baseProc);
}
/*********/
void SFDUtils::
computeIntegratedErrorL2Norm(Real&                                 a_norm,
	                     const Vector<ChomboSFDInterface>&     a_SFDInterface,
	                     const Vector<Real>&                   a_dxVect,
	                     const Vector<int>&                    a_refRatio,
	                     int                                   a_filterIndex,
	                     bool                                  a_incOverlapData)
{

  int finestLevel = a_SFDInterface.size() - 1;

  a_norm = 0.;
  Real locNorm = 0.; // norm on this processor

  Real domainVolume = 0.;
  Real locDomVolume = 0.;

  // First, copy the finest level's data
  // scoping begin
  {   
    Real cellVol = a_dxVect[finestLevel];
    for (int idir = 1; idir < SpaceDim; idir++)
    {
      cellVol *= a_dxVect[finestLevel];
    }

    const LevelData<EBCellFAB>& finestError = *(a_SFDInterface[finestLevel].getIntegratedError(a_filterIndex));

    int nvar = finestError.nComp();
    const DisjointBoxLayout& finestGrids = finestError.getBoxes();
    DataIterator finestDit = finestGrids.dataIterator();
    // iterator over the grids on this processor
    for (finestDit.begin(); finestDit.ok(); ++finestDit)
    {
      const EBCellFAB& error = finestError[finestDit()];
      const Box& box = finestGrids.get(finestDit());
      const EBISBox& ebisBox = error.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        const VolIndex& vof = vofit();
        Real volfrac = ebisBox.volFrac(vof);
        locDomVolume += volfrac*cellVol;

        Real ke = 0.;
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          ke += (error(vofit(), ivar))*(error(vofit(), ivar));
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

    const LevelData<EBCellFAB>& levelError = *(a_SFDInterface[ilev].getIntegratedError(a_filterIndex));

    int nvar = levelError.nComp();
    const DisjointBoxLayout& levelGrids = levelError.getBoxes();
    const DisjointBoxLayout& finerGrids = a_SFDInterface[ilev+1].getIntegratedError(a_filterIndex)->getBoxes(); 

    DataIterator levelDit = levelGrids.dataIterator();
    LayoutIterator finerLit = finerGrids.layoutIterator();

    // iterator over the grids on this processor
    for (levelDit.begin(); levelDit.ok(); ++levelDit)
    {
      const EBCellFAB& error = levelError[levelDit()];
      const Box& thisBox = levelGrids.get(levelDit());
      const EBISBox& ebisBox = error.getEBISBox();
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
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          ke += (error(vofit(), ivar))*(error(vofit(), ivar));
        }
        locNorm += volfrac*cellVol*cellVol*ke;
      }
    }
  }

  int baseProc = 0;
  Vector<Real> volVec;
  Vector<Real> normVec;
  gather(normVec, locNorm, baseProc);
  gather(volVec, locDomVolume, baseProc);

  if (procID() == baseProc)
  {
    CH_assert(normVec.size() == numProc());
    for (int ivec = 0; ivec < normVec.size(); ivec++)
    {
      a_norm += normVec[ivec];
      domainVolume += volVec[ivec];
    }
  }
  //broadcast the sum to all processors.
  broadcast(a_norm, baseProc);
  broadcast(domainVolume, baseProc);

  a_norm = sqrt(a_norm);
  a_norm /= domainVolume;
}
/*********/
