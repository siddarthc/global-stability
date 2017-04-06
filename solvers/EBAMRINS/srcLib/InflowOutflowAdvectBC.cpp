#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBArith.H"
#include "InflowOutflowAdvectBC.H"
#include "DirichletPoissonEBBC.H"

#include <random>

#include "NamespaceHeader.H"

extern Real g_simulationTime;

/****************************/
void
InflowOutflowAdvectBC::
fluxBC(EBFluxFAB&            a_primGdnv,
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
  Box FBox = a_primGdnv[a_dir].getRegion();
  Box cellBox = a_faceBox;
  CH_assert(a_primGdnv[a_dir].nComp() == 1);

  // Determine which side and thus shifting directions
  int isign = sign(a_side);
  cellBox.shiftHalf(a_dir,isign);

  //figure out if we are on a wall and if its an inflow and the inflow value
  bool isInflow   = ((a_side == Side::Lo) && (a_dir==m_flowDir));
  bool isOutflow  = ((a_side == Side::Hi) && (a_dir==m_flowDir));

/*
  //just flipping inflow to outflow in case of adjoint solver
  if (m_adjointSolver)
  {
    bool tmp = isInflow;
    isInflow = isOutflow;
    isOutflow = tmp;
  }
*/

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(cellBox))
  {
    cellBox &= m_domain;
    // Find the strip of cells next to the domain boundary
    Box bndryBox = adjCellBox(cellBox, a_dir, a_side, 1);

    // Shift things to all line up correctly
    bndryBox.shift(a_dir,-isign);

    IntVectSet ivs(bndryBox);
    for (VoFIterator vofit(ivs, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Vector<FaceIndex> bndryFaces = a_ebisBox.getFaces(vof, a_dir, a_side);
      for (int iface= 0; iface < bndryFaces.size(); iface++)
      {
        //this all will work for the case of spatially varying inflow
        //for inhomogeneous other faces, this is wrong.
        const FaceIndex& face = bndryFaces[iface];
        Real velComp = 1.e99;

        if (!m_adjointSolver)
        {
          if (isInflow && (m_velComp == m_flowDir))
          {
            if (!m_doPoiseInflow && !m_doWomersleyInflow)
            {
              velComp = m_inflowVel;
            }
            else if (m_doPoiseInflow)
            {
              //set coordinates based on slip walls
              RealVect prob_lo = RealVect::Zero;
              //assumes inflow boundary is always on lo side of domain
              const RealVect loc  = EBArith::getFaceLocation(face, m_dx, prob_lo);
              Real radius = m_poiseInflowFunc->getRadius(loc);
              velComp = m_poiseInflowFunc->getVel(radius)[m_flowDir];
            }
            else if (m_doWomersleyInflow)
            {
              double PI = 3.1416;
              int freq[10] =
              {
                1,2,3,4,5,6,7,8
              };
              double Vp[10] =
              {
                0.33,0.24,0.24,0.12,0.11,0.13,0.06,0.04
              };
              double Theta[10] =
              {
                74,79,121,146,147,179,233,218
              };
              double AmpXsi[10] =
              {
                1.7639,1.4363,1.2517,1.1856,1.1603,1.1484,1.1417,1.1214
              };
              double AngXsi[10] =
              {
                -0.2602,-0.3271,-0.2799,-0.2244,-0.1843,-0.1576,-0.1439,-0.1195
              };
  
              Real VelMult = 2;

              int Maxp = 8;

              for (int p=0;p<Maxp;p++)
                VelMult += Vp[p] * AmpXsi[p] * cos (2 * PI * freq[p] * g_simulationTime - Theta[p]*PI/180 + AngXsi[p]);

              velComp = m_inflowVel*VelMult/2;
            }
          }
          else if (isInflow && (m_velComp != m_flowDir))
          {
            velComp = 0.0;
          }
          else if (isInflow && m_params.m_doInflowWhiteNoise)
          {
            std::default_random_engine generator;
            std::normal_distribution<double> distribution(0.0,m_params.m_inflowSigma);
            velComp += distribution(generator);
          }
          else if (isOutflow)
          {
            velComp = a_primExtrap(vof, 0);
          }
          else //solid wall
          {
            if (a_dir == m_velComp)
            {
              velComp = 0.0;
            }
            else
            {
              velComp = a_primExtrap(vof, 0);
            }
          }
        } // end !adjointSolver
        else  // adjointSolver
        {
          if (!isOutflow)
          {
            velComp = 0.; // homogeneous Dirichlet
          }
          else
          {
            // set Ubase.n*vel + 1/Re * dvel/dn = 0
            CH_assert(m_RobinBCData != NULL);
            // dvel/dn = primExtrap - primCenter
            int lev = (*m_RobinBCData).levelDxMap.find(m_dx[0])->second;
            Real baseDataVal = (*((*(*m_RobinBCData).dataPtr)[lev]))[a_dit][a_dir](face, 0); 
            baseDataVal *= -1.*(*m_RobinBCData).Re;
            Real factor = baseDataVal/(baseDataVal - isign*8./(3.*m_dx[a_dir]));
//            velComp = (a_primExtrap(vof, 0) - a_primCenter(vof, 0))/(-1.*baseDataVal*(*m_RobinBCData).Re*0.5*m_dx[0]);

            RealVect prob_lo = RealVect::Zero;
            RobinPoissonBCFunctions::getHigherOrderFaceFlux(velComp, vof, 0, a_primCenter, prob_lo, m_dx, a_dir, a_side, a_dit, a_time, false, *m_RobinBCData);
            // velComp comtained dvel/dn = C*velOutflow
            velComp /= factor;
          }
        }
        a_primGdnv[a_dir](face, 0) = velComp;
      }
    }
  }
}

#include "NamespaceFooter.H"
