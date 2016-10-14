/*
 * This pig is added to protect this code from the bugs.
 *      To make best use of this pig, compile with opt and fPIG flags
 *                               _
 *       _._ _..._ .-',     _.._(`))
 *      '-. `     '  /-._.-'    ',/
 *         )         \            '.
 *        / _    _    |             \
 *       |  o    o    /              |
 *       \   .-.                     ;
 *        '-('' ).-'       ,'       ;
 *           '-;           |      .'
 *              \           \    /
 *              | 7  .__  _.-\   \
 *              | |  |  ``/  /`  /
 *             /,_|  |   /,_/   /
 *                /,_/      '`-'
 *
 *
 * refer = curious ?
 *         https://gist.github.com/crittelmeyer/b4c8ec6f02c024798ccb
 *         : NULL;
**/

#include <cmath>

#include "InflowOutflowIBC.H"
#include "EBISLayout.H"
#include "EBLoHiCenter.H"
#include "ParmParse.H"
#include "VoFIterator.H"
#include "EBLGIntegrator.H"
#include "EBLevelDataOps.H"
#include "EBPatchPolytropicF_F.H"
#include "EBSolidF_F.H"
#include "BaseFab.H"

#include <random>

#include "NamespaceHeader.H"

/****************************/
/****************************/
InflowOutflowIBC::
InflowOutflowIBC(const InflowOutflowParams&  a_params)
{
  m_params = a_params;
  FORT_SETGAMMAANDSMALL(CHF_CONST_REAL(a_params.m_gamma));

  m_isDefined = false;
}
/****************************/
/****************************/
void
InflowOutflowIBC::
setBndrySlopes(EBCellFAB&       a_deltaPrim,
               const EBCellFAB& a_primState,
               const EBISBox&   a_ebisBox,
               const Box&       a_box,
               const int&       a_dir)
{
}

/****************************/
/****************************/
void
InflowOutflowIBC::
fluxBC(EBFluxFAB&            a_flux,
       const EBCellFAB&      a_primCenter,
       const EBCellFAB&      a_primExtrap,
       const Side::LoHiSide& a_side,
       const Real&           a_time,
       const EBISBox&        a_ebisBox,
       const DataIndex&      a_dit,
       const Box&            a_box,
       const Box&            a_faceBox,
       const int&            a_dir)
{
  CH_assert(m_isDefined);
  CH_assert(!m_domain.isPeriodic(a_dir));

  //gets face centered region of data
  Box FBox = a_flux[a_dir].getSingleValuedFAB().box();
  Box cellBox = FBox;
  int numFlux = a_flux[a_dir].nComp();

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  //figure out if we are on a wall and if its an inflow and the inflow value
  bool isInflow   = ((a_side == Side::Lo) && (a_dir==m_params.m_inflowDir));
  bool isOutflow  = ((a_side == Side::Hi) && (a_dir==m_params.m_outflowDir));
  bool isExtrap   = ((!(isInflow)) && (!(isOutflow)) && (((a_side == Side::Lo) && (m_params.m_wallBCLo[a_dir]==0)) || ((a_side == Side::Hi) && (m_params.m_wallBCHi[a_dir]==0))));
  // else hardwiring it to be wall

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
    {
      cellBox &= m_domain;
      // Find the strip of cells next to the domain boundary
      Box boundaryBox = bdryBox(cellBox, a_dir, a_side, 1);

      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-isign);
      IntVectSet ivs(boundaryBox);
      VoFIterator vofit(ivs, a_ebisBox.getEBGraph());
      // Set the boundary fluxes
      for(vofit.reset(); vofit.ok(); ++vofit)
        {
          const VolIndex& vof = vofit();
          Vector<Real> qgdnv(QNUM);
          Vector<Real> fluxv(FNUM);
          Vector<Real> qgdnv1(QNUM);
          Vector<Real> qgdnv2(QNUM);

          if(isInflow)  //inflow side
            {
              int inormVelVar = CMOMX + a_dir;
              Real u_o = isign*m_params.m_inflowVel[a_dir];
              Real p_o = m_params.m_inflowPress;
              Real rho_o = m_params.m_inflowDense;
              Real T_o = m_params.m_inflowTemp;
              Real c_o = sqrt(m_params.m_gamma*p_o/rho_o);

              std::default_random_engine generator;
              std::normal_distribution<double> distribution(0.0,m_params.m_inflowSigma);

              Real u_i = isign*a_primExtrap(vof, inormVelVar);
              Real p_i = a_primExtrap(vof, QPRES);
              Real rho_i = a_primExtrap(vof, QRHO);
              Real T_i = a_primExtrap(vof, QCVTEMP);
              Real c_i = sqrt(m_params.m_gamma*p_i/rho_i);

              if (abs(u_o) < c_o)
                {
                  if (m_params.m_verbosity >= 3)
                    {
                      pout() << "subsonic inlet : ";
                    }

                  if (m_params.m_totPresBC)
                    {
                      if (m_params.m_verbosity >= 3)
                       {
                         pout() << "Using total pressure BC..." << endl;
                       }
                      Real Rplus = u_i + 2.*c_i/(m_params.m_gamma -1.);
                      Real Ht = c_i*c_i/(m_params.m_gamma - 1.) + 0.5*(Rplus + 2.*c_i/(m_params.m_gamma-1.))*(Rplus + 2.*c_i/(m_params.m_gamma-1.));

                      Real a = 1. + 2./(m_params.m_gamma-1.);
                      Real b = 2.*Rplus;
                      Real c = 0.5*(m_params.m_gamma-1.)*(Rplus*Rplus - 2.*Ht);

                      Real cb1 = -0.5*b/a + 0.5*(sqrt(b*b - 4.*a*c))/a;
                      Real cb2 = -0.5*b/a - 0.5*(sqrt(b*b - 4.*a*c))/a;

                      Real cb = max(cb1, cb2);

                      Real ub = 2.*cb/(m_params.m_gamma-1.)-Rplus;
                      Real Mb = ub/cb;

                      Real M_set = u_o/c_o;
                      Real Pt_set = p_o*pow((1. + 0.5*M_set*M_set*(m_params.m_gamma -1.)) , (m_params.m_gamma/(m_params.m_gamma-1.)));
                      Real T_o = m_params.m_inflowTemp;
                      Real Tt_set = T_o*(1. + M_set*M_set*0.5*(m_params.m_gamma - 1.));
                      Real Pb = Pt_set/pow((1. + 0.5*Mb*Mb*(m_params.m_gamma -1.)) , (m_params.m_gamma/(m_params.m_gamma-1.)));
                      Real Tb = Tt_set/(1. + Mb*Mb*0.5*(m_params.m_gamma - 1.));
                      Real Rhob = Pb/(Tb*(m_params.m_specificHeat)*(m_params.m_gamma-1.));
                      qgdnv[QRHO] = Rhob;
                      qgdnv[QPRES] = Pb;

                      for(int idir = 0; idir < SpaceDim; idir++)
                        {
                          qgdnv[QVELX+idir] = m_params.m_inflowVel[idir] + distribution(generator);
                          if (idir == a_dir) qgdnv[QVELX+idir] = ub + distribution(generator);
                        }
                    }
                  else
                    {
                      if (m_params.m_verbosity >= 3)
                        {
                          pout() << "Using Toro BC..." << endl;
                        }
                      for(int idir = 0; idir < SpaceDim; idir++)
                        {
                          qgdnv[QVELX+idir] = m_params.m_inflowVel[idir] + distribution(generator);
                          qgdnv1[QVELX+idir] = a_primExtrap(vof, QVELX+idir);
                        }
                      qgdnv[QPRES] = m_params.m_inflowPress;
                      qgdnv1[QPRES] = a_primExtrap(vof, QPRES);
                      qgdnv[QRHO] = m_params.m_inflowDense;
                      qgdnv1[QRHO] = a_primExtrap(vof, QRHO);
                      qgdnv[QCVTEMP] = m_params.m_inflowTemp*m_params.m_specificHeat;
                      qgdnv1[QCVTEMP] = a_primExtrap(vof, QCVTEMP);
                    }
                }
              else
                {
                  if (m_params.m_verbosity >= 3)
                    {
                      pout() << "supersonic inlet" << endl;
                    }
                  qgdnv[QRHO] = m_params.m_inflowDense;
                  for(int idir = 0; idir < SpaceDim; idir++)
                    {
                      qgdnv[QVELX+idir] = m_params.m_inflowVel[idir] + distribution(generator);
                    }
                  qgdnv[QPRES] = m_params.m_inflowPress;
                  qgdnv[QCVTEMP] = 0.; // just a dummy value
                }
            }
          else if (isOutflow || isExtrap)//outflow side --extrapolation
            {
              for(int ivar = 0; ivar < QNUM; ivar++)
                {
                  qgdnv[ivar] = a_primExtrap(vof, ivar);
                }
            }
        Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
        for(int iface= 0; iface < bndryFaces.size(); iface++)
          {
            const FaceIndex& face = bndryFaces[iface];
            if (isInflow)
              {
                Real u = qgdnv[QVELX];
                Real p = qgdnv[QPRES];
                Real rho = qgdnv[QRHO];
                Real c = m_params.m_gamma*p/rho;

                //solve riemann problem at inlet only for subsonic flow
                if ((abs(u) < c) && !(m_params.m_totPresBC) )
                  {
                    /**/
                    FORT_POINTRIEMANN(CHF_VR(qgdnv),
                                      CHF_VR(qgdnv1),
                                      CHF_VR(qgdnv2),
                                      CHF_CONST_INT(a_dir));
                    /**/
                    FORT_POINTGETFLUX(CHF_VR(fluxv),
                                      CHF_VR(qgdnv2),
                                      CHF_CONST_INT(a_dir));
// NOTE:- need to figure out a fix to make this work for flow coming in 
//        from the high side.
                  }
                else
                  {
                    /**/
                    FORT_POINTGETFLUX(CHF_VR(fluxv),
                                      CHF_VR(qgdnv),
                                      CHF_CONST_INT(a_dir));
                  }
                for(int ivar = 0; ivar < FNUM; ivar++)
                  {
                    a_flux[a_dir](face, ivar) = fluxv[ivar];
                  }
              }
            else if (isOutflow || isExtrap)
              {
                Real u = qgdnv[QVELX];
                Real p = qgdnv[QPRES];
                Real rho = qgdnv[QRHO];
                Real c = m_params.m_gamma*p/rho;
                if (u < c)
                  {
                    /**/
                    FORT_POINTRIEMANN(CHF_VR(qgdnv),
                                      CHF_VR(qgdnv),
                                      CHF_VR(qgdnv1),
                                      CHF_CONST_INT(a_dir));
                    /**/
                    FORT_POINTGETFLUX(CHF_VR(fluxv),
                                      CHF_VR(qgdnv1),
                                      CHF_CONST_INT(a_dir));
                    /**/
                  }
                else
                  {
                    FORT_POINTGETFLUX(CHF_VR(fluxv),
                                      CHF_VR(qgdnv),
                                      CHF_CONST_INT(a_dir));
                  }
                for (int ivar = 0; ivar < FNUM; ivar++)
                  {
                    a_flux[a_dir](face, ivar) = fluxv[ivar];
                  }
              }
            else // hardwiring to be a wall
              {
                //set all fluxes to zero then fix normal momentum
                for(int ivar = 0; ivar < numFlux; ivar++)
                  {
                    a_flux[a_dir](face, ivar) = 0;
                  }
                int inormVelVar = CMOMX + a_dir;
                int inormMomVar = CMOMX + a_dir;
                Real press = a_primExtrap(vof, QPRES);
                Real dense = a_primExtrap(vof, QRHO);
                Real unorm = a_primExtrap(vof, inormVelVar);
                Real speed = sqrt(m_params.m_gamma*press/dense);


                a_flux[a_dir](face, inormMomVar) = press + isign*dense*unorm*speed;
              }
          }
      }
  }
}
/****************************/
/****************************/
void
InflowOutflowIBC::define(const ProblemDomain&  a_domain,
                       const RealVect& a_dx)
{
  m_domain = a_domain;
  m_dx = a_dx[0];
  m_isDefined = true;
}
/****************************/
/****************************/
void
InflowOutflowIBC::
initialize(LevelData<EBCellFAB>& a_conState,
           const EBISLayout& a_ebisl) const
{
  CH_assert(m_isDefined);
  Real initialR = m_params.m_initDense;
  RealVect initialMom = m_params.m_initVel;
  initialMom *= initialR;
  Real initKE = 0.;
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      initKE += 0.5*(m_params.m_initVel[idir])*(m_params.m_initVel[idir]);
    }

  initKE *= initialR;
  Real initIE = initialR*m_params.m_specificHeat*m_params.m_initTemp;

  Real initialE = initKE + initIE;

  if (!m_params.m_doInitWhiteNoise)
    {
      for(int idir = 0; idir < SpaceDim; idir++)
        {
          EBLevelDataOps::setVal(a_conState, initialMom[idir], CMOMX+idir);
        }
      EBLevelDataOps::setVal(a_conState, initialR, CRHO);
      EBLevelDataOps::setVal(a_conState, initialE, CENG);
    }
  else
    {
      std::default_random_engine generator;
      std::normal_distribution<double> distribution(0.0,m_params.m_initSigma);
      EBLevelDataOps::setVal(a_conState, initialR, CRHO);
      EBLevelDataOps::setVal(a_conState, initialE, CENG);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
            EBLevelDataOps::setVal(a_conState, initialMom[idir], CMOMX+idir);
        }

      for (DataIterator dit = a_conState.dataIterator(); dit.ok(); ++dit)
        {
          DataIndex d = dit();
          EBCellFAB& conState = a_conState[d];
          BaseFab<Real>& regCons = conState.getSingleValuedFAB();
          const EBISBox& ebisBox = conState.getEBISBox();
          const Box& region = conState.getRegion();
          for (BoxIterator bit(region); bit.ok(); ++bit)
            {
              const IntVect& iv = bit();
              if (ebisBox.isRegular(iv))
                {
                  for (int idir = 0; idir < SpaceDim; idir++)
                    {
                      Real randVal = distribution(generator);
                      regCons(iv,CMOMX+idir) += randVal;
                    }
                }
            }
        }
    }
}
/****************************/
/****************************/
InflowOutflowIBC::
~InflowOutflowIBC()
{
}
/****************************/
/****************************/

#include "NamespaceFooter.H"
