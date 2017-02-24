/*
 *
 *
 *
 *
 */

#include "EBAMRLinINS.H"
#include "ParmParse.H"
#include "EBAMRNoSubcycleF_F.H"
#include "EBAMRPoissonOpF_F.H"
#include "MeshRefine.H"
#include "BRMeshRefine.H"
#include "EBEllipticLoadBalance.H"
#include "EBArith.H"
#include "EBPWLFineInterp.H"
#include "EBCoarseAverage.H"
#include "EBFluxFactory.H"
#include "EBCellFactory.H"
#include "EBLevelAdvect.H"
#include "EBGradDivFilter.H"
#include "EBPatchAdvect.H"
#include "REAL.H"
#include "EBPhysIBCFactory.H"
#include "EBAMRIO.H"
#include "BaseIFFactory.H"
#include "EBLevelRedist.H"
#include "BaseIVFactory.H"
#include "EBConstantCFInterp.H"
#include "EBArith.H"
#include "EBAMRDataOps.H"
#include "NeumannPoissonEBBC.H"
#include "DirichletPoissonEBBC.H"
#include "InflowOutflowIBC.H"
#include "EBNormalizeByVolumeFraction.H"
#include <iomanip>
#include <cmath>
#include <cstdio>
#include "memusage.H"
#include "memtrack.H"

extern Real g_simulationTime;
#define debugIV IntVect(D_DECL(16, 3, 0))

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/
EBAMRLinINS::
EBAMRLinINS(const AMRParameters&      a_params,
            const EBIBCFactory&       a_ibcfact,
            const ProblemDomain&      a_coarsestDomain,
            Real                      a_viscosity,
            bool                      a_doAdjoint,
            const EBIndexSpace* const a_ebisPtr):
  m_ebisPtr(a_ebisPtr)
{
  if (a_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::EBAMRLinINS" << endl;
    }

  //set parameters of the run
  m_params    = a_params;
  m_viscosity = a_viscosity;
  m_viscousCalc = (m_viscosity > 0);

  //create initial and boundary condition object
  m_ibc    =   a_ibcfact.create();

  //resize vectors and set them where we can
  Real coarsestDx = m_params.m_domainLength/Real(a_coarsestDomain.size(0));
  int nlevels = m_params.m_maxLevel + 1;
  m_domain.resize(nlevels);
  m_dx.resize(nlevels);
  m_grids.resize(nlevels);

  m_ebisl.resize(nlevels);
  m_eblg.resize(nlevels);

  m_quadCFI.resize(nlevels);
  m_aveOper.resize(nlevels);
  m_aveSpac.resize(nlevels);
  m_ebLevAd.resize(nlevels);
  m_fluxReg.resize(nlevels);
  m_velo.resize(nlevels, NULL);
  m_pres.resize(nlevels, NULL);
  m_gphi.resize(nlevels, NULL);

  m_coveredFaceLitLo.resize(nlevels, NULL);
  m_coveredFaceLitHi.resize(nlevels, NULL);
  m_coveredSetsLitLo.resize(nlevels, NULL);
  m_coveredSetsLitHi.resize(nlevels, NULL);

  allocateDataHolders();

  m_domain[0] = a_coarsestDomain;
  m_dx[0]     =   coarsestDx;
  for (int ilev = 1; ilev < nlevels; ilev++)
    {
      CH_assert(m_params.m_refRatio[ilev-1] > 0);
      m_domain[ilev] = refine(m_domain[ilev-1], m_params.m_refRatio[ilev-1]);
      m_dx[ilev] = m_dx[ilev-1]/Real(m_params.m_refRatio[ilev-1]);
    }
  m_prescribedDt = -1.0;
  m_useFixedDt = false;
  m_steadyState = false;
  m_stopAdvance = false;

  m_ccProjector  = NULL;
  m_macProjector = NULL;
  m_time = 0.0;
  m_curStep = 0;
  m_dt = -1.0;

  //setup still needs to get called
  m_isSetup  = false;
  m_pointsUpdated = 0;
}
/*********/
void EBAMRLinINS::
allocateDataHolders()
{
  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      m_velo[ilev]   = new LevelData<EBCellFAB>();
      m_pres[ilev]   = new LevelData<EBCellFAB>();
      m_gphi[ilev]   = new LevelData<EBCellFAB>();
 
      m_coveredFaceLitLo[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredFaceLitHi[ilev] = new LayoutData< Vector< Vector<VolIndex> > >();
      m_coveredSetsLitLo[ilev] = new LayoutData< Vector< IntVectSet > >      ();
      m_coveredSetsLitHi[ilev] = new LayoutData< Vector< IntVectSet > >      ();
    }
}
/*********/
EBAMRLinINS::
~EBAMRLinINS()
{
  if (m_params.m_verbosity > 3)
    {
      pout() << "EBAMRLinINS::~EBAMRLinINS" << endl;
    }
  delete m_ibc;

  for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
      delete m_velo[ilev];
      delete m_pres[ilev];
      delete m_gphi[ilev];

      delete m_coveredFaceLitLo[ilev];
      delete m_coveredFaceLitHi[ilev];
      delete m_coveredSetsLitLo[ilev];
      delete m_coveredSetsLitHi[ilev];
  
    }

  if (m_ccProjector !=  NULL)
    {
      delete m_ccProjector;
    }
  if (m_macProjector !=  NULL)
    {
      delete m_macProjector;
    }
}
