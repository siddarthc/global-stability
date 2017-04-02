/*
*
*
*
*/

#include "WavemakerUtils.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"
#include "EBLevelDataOps.H"
#include "EBAMRDataOps.H"
#include "EBAMRIO.H"
#include <cmath>

/*********/
void WavemakerUtils::
plotWavemaker(std::string               a_fileName,
              std::string               a_directModeFile,
              std::string               a_adjointModeFile,
              Vector<int>               a_refRatio,
              const ProblemDomain&      a_coarsestDomain,
              Real                      a_coarsestDx,
              const EBIndexSpace* const a_ebisPtr) 
{
  Vector<LevelData<EBCellFAB>* > directModeVelocity;
  Vector<LevelData<EBCellFAB>* > adjointModeVelocity;
  Vector<LevelData<EBCellFAB>* > structuralSensitivity;

  Vector<DisjointBoxLayout> grids;
  Vector<EBLevelGrid> eblgs;

  int nlevels;

#ifdef CH_USE_HDF5

  { // scoping
    HDF5Handle handleIn(a_directModeFile, HDF5Handle::OPEN_RDONLY);
    HDF5HeaderData header;
    header.readFromFile(handleIn);

    int finestLevel = header.m_int["finest_level"];
    nlevels = finestLevel+1;

    directModeVelocity.resize(nlevels);
    adjointModeVelocity.resize(nlevels);
    structuralSensitivity.resize(nlevels);
    grids.resize(nlevels);
    eblgs.resize(nlevels);

    int nGhost = 1.;
    int nEBGhost = 1.;
    int nComp = SpaceDim;

    // get all the grids 
    for (int ilev=0; ilev <= finestLevel; ilev++)
    {
      directModeVelocity[ilev] = new LevelData<EBCellFAB>();
      adjointModeVelocity[ilev] = new LevelData<EBCellFAB>();
      structuralSensitivity[ilev] = new LevelData<EBCellFAB>();

      handleIn.setGroupToLevel(ilev);
      // Get the grids
      Vector<Box> vboxGrids;
      const int gridStatus = read(handleIn, vboxGrids);
      if (gridStatus != 0)
      {
        MayDay::Error("readCheckpointLevel: file has no grids");
      }

      Vector<int> proc_map;
      ProblemDomain levelDomain = a_coarsestDomain;
      if (ilev > 0)
      {
        CH_assert(a_refRatio.size() >= ilev-1);
        for (int ilev1 = 1; ilev1 <= ilev; ilev1++)
        {
          levelDomain.refine(a_refRatio[ilev1-1]);
        }
      }

      LoadBalance(proc_map, vboxGrids);
      broadcast(proc_map, uniqueProc(SerialTask::compute));
      grids[ilev] = DisjointBoxLayout(vboxGrids, proc_map);
      eblgs[ilev] = EBLevelGrid(grids[ilev], levelDomain, nEBGhost, a_ebisPtr);
    
      EBCellFactory ebcellfact(eblgs[ilev].getEBISL());
      directModeVelocity[ilev]->define(grids[ilev], nComp, IntVect::Unit, ebcellfact);
      adjointModeVelocity[ilev]->define(grids[ilev], nComp, IntVect::Unit, ebcellfact);
      structuralSensitivity[ilev]->define(grids[ilev], 1, IntVect::Unit, ebcellfact);

      EBLevelDataOps::setToZero(*directModeVelocity[ilev]);
      EBLevelDataOps::setToZero(*adjointModeVelocity[ilev]);
      EBLevelDataOps::setToZero(*structuralSensitivity[ilev]);
    
      read<EBCellFAB>(handleIn, *directModeVelocity[ilev], "velo", grids[ilev], Interval(), false);
    }

    handleIn.close();
  } // end scoping
  
  { //scoping
    HDF5Handle handleIn(a_adjointModeFile, HDF5Handle::OPEN_RDONLY);  
    HDF5HeaderData header;
    header.readFromFile(handleIn);

    int finestLevel = header.m_int["finest_level"];
    CH_assert(finestLevel+1 == nlevels);
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
      handleIn.setGroupToLevel(ilev);
      read<EBCellFAB>(handleIn, *adjointModeVelocity[ilev], "velo", grids[ilev], Interval(), false);
    }

    handleIn.close();
  } // end scoping

  // compute structural sensitivity
  computeStructuralSensitivity(structuralSensitivity, 
                               directModeVelocity, 
                               adjointModeVelocity,
                               a_refRatio,
                               a_coarsestDx);

  EBAMRDataOps::averageDown(structuralSensitivity, eblgs, a_refRatio);

  // set covered stuff to zero
  EBAMRDataOps::setCoveredVal(structuralSensitivity, 0.);

  { // plot structural sensitivity
    Vector<Real> coveredValues;
    Vector<string> varNames(1, string("structuralSensitivity"));
    std::string name = "plot_" + a_fileName;
    writeEBHDF5(name, 
                grids, 
                structuralSensitivity, 
                varNames, 
                a_coarsestDomain.domainBox(),
                a_coarsestDx,
                0.,
                0.,
                a_refRatio,
                nlevels,
                false,
                coveredValues);
  }

  { // checkpoint structural sensitivity
    HDF5HeaderData header;
    header.m_int["finest_level"] = nlevels-1;

    std::string name = "check_" + a_fileName;
    HDF5Handle handleOut(name, HDF5Handle::CREATE);
    //Write the header for this level
    header.writeToFile(handleOut);
    for (int ilev = 0; ilev < nlevels; ilev++)
    {
      handleOut.setGroupToLevel(ilev);
      write(handleOut,grids[ilev]);
      write(handleOut,*structuralSensitivity[ilev],"structuralSensitivity");
    }

    handleOut.close();
  }

  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    delete directModeVelocity[ilev];
    delete adjointModeVelocity[ilev];
    delete structuralSensitivity[ilev];
  }

#else

  MayDay::Error("Chombo needs HDF5 to read baseflowFile");

#endif
}
/*********/
void WavemakerUtils::
computeStructuralSensitivity(Vector<LevelData<EBCellFAB>* >& a_out,
                             const Vector<LevelData<EBCellFAB>* >& a_directModeVelo,
                             const Vector<LevelData<EBCellFAB>* >& a_adjointModeVelo,
                             Vector<int> a_refRatio,
                             Real a_coarsestDx)
{
  int nlevels = a_out.size();
  Real dotProd = 0.;
  computeKappaDotProduct(dotProd, a_directModeVelo, a_adjointModeVelo, a_refRatio, a_coarsestDx);  

  for (int ilev = 0; ilev < nlevels; ilev++)
  {
    LevelData<EBCellFAB>& levelData = *a_out[ilev];
    const DisjointBoxLayout& levelGrids = levelData.getBoxes();
    DataIterator dit = levelGrids.dataIterator();
    // iterate over grids on this processor
    for (dit.begin(); dit.ok(); ++dit)
    {
      EBCellFAB& outFAB = levelData[dit()];
      const EBCellFAB& dirVelFAB = (*a_directModeVelo[ilev])[dit()];
      const EBCellFAB& adjVelFAB = (*a_adjointModeVelo[ilev])[dit()];
      const Box& box = levelGrids.get(dit());
      const EBISBox& ebisBox = outFAB.getEBISBox();
      IntVectSet ivs(box);
      for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
      {
        Real ke1 = 0.;
        Real ke2 = 0.;
        for (int idir = 0; idir < SpaceDim; idir++)
        {
          ke1 += dirVelFAB(vofit(), idir)*dirVelFAB(vofit(), idir);
          ke2 += adjVelFAB(vofit(), idir)*adjVelFAB(vofit(), idir);
        }
        outFAB(vofit(), 0) = sqrt(ke1)*sqrt(ke2);
        outFAB(vofit(), 0) /= abs(dotProd);
      }
    }
  }

  Real max, min;
  EBAMRDataOps::getMaxMin(max, min, a_out, 0);
  EBAMRDataOps::scale(a_out, 1./max);
}
/*********/
void WavemakerUtils::
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
