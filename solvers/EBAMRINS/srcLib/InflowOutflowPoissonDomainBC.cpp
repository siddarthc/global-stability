#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "InflowOutflowPoissonDomainBC.H"
#include "PoiseuilleInflowBCValue.H"
#include "NeumannPoissonDomainBC.H"
#include "DirichletPoissonDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletPoissonEBBC.H"
#include "VoFIterator.H"

#include "NamespaceHeader.H"

extern Real g_simulationTime;
bool s_higherOrderHelmBC = false;

//pressure stencil
void InflowOutflowPoissonDomainBC::getFluxStencil(      VoFStencil&      a_stencil,
                                                  const VolIndex&        a_vof,
                                                  const int&             a_comp,
                                                  const RealVect&        a_dx,
                                                  const int&             a_idir,
                                                  const Side::LoHiSide&  a_side,
                                                  const EBISBox&         a_ebisBox)
{
  bool isOutflow= ((a_side==Side::Hi) && (a_idir==m_flowDir));
  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.setEBOrder(1);
      diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
}
/*********/
//called by MAC projector enforceFaceVel
void
InflowOutflowPoissonDomainBC::
enforceFaceVel(LevelData<EBFluxFAB>&    a_velocity,
               const DisjointBoxLayout& a_grids,
               const EBISLayout&        a_ebisl,
               const ProblemDomain&     a_domain,
               const RealVect&          a_dx,
               const Real&              a_time,
               const RealVect&          a_origin,
               const bool&              a_doDivFreeOutflow)
{
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          DirichletPoissonEBBC::s_velComp = idir;
          if (!a_domain.isPeriodic(idir))
            {
              EBFluxFAB& fluxVel = a_velocity[dit()];
              CH_assert(fluxVel.nComp() == 1);
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Box sideBox = adjCellBox(grid, idir, sit(), 1);
                  int ishift = -sign(sit());
                  sideBox.shift(idir, ishift);
                  IntVectSet ivsIrreg(sideBox);
                  FaceStop::WhichFaces stopCrit;
                  if (sit() == Side::Lo)
                    {
                      stopCrit = FaceStop::LoBoundaryOnly;
                    }
                  else
                    {
                      stopCrit = FaceStop::HiBoundaryOnly;
                    }

                  FaceIterator faceit(ivsIrreg ,ebgraph, idir, stopCrit);
                  for (faceit.reset(); faceit.ok(); ++faceit)
                    {
                      Real boundaryVel = 0.0;
                      getFaceVel(boundaryVel, faceit(), fluxVel, a_origin, a_dx, idir, 0, a_time,sit(),a_doDivFreeOutflow, dit());
                      fluxVel[idir](faceit(), 0) = boundaryVel;
                    }
                }
            }
        }
    }
  a_velocity.exchange();
}

/*********/
void
InflowOutflowPoissonDomainBC::
enforceFaceVel(LevelData<EBFluxFAB>&    a_velocity,
               const DisjointBoxLayout& a_grids,
               const EBISLayout&        a_ebisl,
               const ProblemDomain&     a_domain,
               const RealVect&          a_dx,
               const Real&              a_time,
               const RealVect&          a_origin,
               const bool&              a_doDivFreeOutflow,
               const int&               a_comp)
{
  DirichletPoissonEBBC::s_velComp = a_comp;
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& grid = a_grids.get(dit());
      const EBGraph& ebgraph = a_ebisl[dit()].getEBGraph();
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (!a_domain.isPeriodic(idir))
            {
              EBFluxFAB& fluxVel = a_velocity[dit()];
              for (SideIterator sit; sit.ok(); ++sit)
                {
                  Box sideBox = adjCellBox(grid, idir, sit(), 1);
                  int ishift = -sign(sit());
                  sideBox.shift(idir, ishift);
                  IntVectSet ivsIrreg(sideBox);
                  FaceStop::WhichFaces stopCrit;
                  if (sit() == Side::Lo)
                    {
                      stopCrit = FaceStop::LoBoundaryOnly;
                    }
                  else
                    {
                      stopCrit = FaceStop::HiBoundaryOnly;
                    }

                  FaceIterator faceit(ivsIrreg ,ebgraph, idir, stopCrit);
                  for (faceit.reset(); faceit.ok(); ++faceit)
                    {
                      Real boundaryVel = 0.0;
                      getFaceVel(boundaryVel, faceit(), fluxVel, a_origin, a_dx, idir, 0, a_time,sit(), a_doDivFreeOutflow, dit());
                      fluxVel[idir](faceit(), 0) = boundaryVel;
                    }
                }
            }
        }
    }
  a_velocity.exchange();
}
/*********/
///
/**
   This is called by cc projection to enforceVelocityBCs on face
   velocities averaged from centers before calculating divergence in mac projection.
*/
void
InflowOutflowPoissonDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  CH_assert(a_vel.nComp() == 1);

  bool isInflow = (a_side==Side::Lo) && (a_idir==m_flowDir);
  bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);

  RealVect point =  EBArith::getFaceLocation(a_face, a_dx, a_probLo);
  RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

  const EBISBox& ebisBox= a_vel[0].getEBISBox();

  if (!m_adjointSolver)
  {

/*
  // flip inflow and outflow in case of adjoint
  if (m_adjointSolver)
  {
    bool tmp = isInflow;
    isInflow = isOutflow;
    isOutflow = tmp;
  }
*/

//    RealVect point =  EBArith::getFaceLocation(a_face, a_dx, a_probLo);
//    RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

//    const EBISBox& ebisBox= a_vel[0].getEBISBox();

    if (isOutflow)
    {
      //quadratic extrapolation with homogeneous Neumann for outflow bc
      a_faceFlux = EBArith::extrapFaceVelToOutflow(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0);//be careful of this 0, it is the component
    }
    else if (isInflow && (velcomp==m_flowDir))
    {
      //input component is always zero.
      //value of face direction corresponds to velocity direction here
      if (!m_doPoiseInflow && !m_doWomersleyInflow)
      {
        a_faceFlux = m_inflowVel;
      }
      else if (m_doPoiseInflow)
      {
        RealVect prob_lo = RealVect::Zero;
        const RealVect loc  = EBArith::getFaceLocation(a_face, a_dx, prob_lo);
        Real radius = m_poiseInflowFunc->getRadius(loc);
        a_faceFlux = m_poiseInflowFunc->getVel(radius)[m_flowDir];
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

        a_faceFlux = m_inflowVel*VelMult/2;
      }
    }
    else
    {
      //must be a solid wall or tangential velocity at inflow
      a_faceFlux = 0;
    }
  } // end !adjointSolver
  else
  {
//    a_faceFlux = RobinBCFunctions::extrapToBoundary(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0,*m_RobinBCData,a_dx[0],a_dit);

    MayDay::Error("InflowOutflowPoissonDomainBC::getFaceVel wrong call for adjoint solver");
  }
}

/*********/
void
InflowOutflowPoissonDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow,
           const DataIndex&      a_dit)
{
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  CH_assert(a_vel.nComp() == 1);

  if (m_adjointSolver)
  {
    bool isOutflow= (a_side==Side::Hi) && (a_idir==m_flowDir);
    const EBISBox& ebisBox= a_vel[0].getEBISBox();
    if (!isOutflow)
    {
      a_faceFlux = 0.;
    }
    else
    {
      a_faceFlux = RobinPoissonBCFunctions::extrapToBoundary(a_face,a_side,a_idir,ebisBox.getEBGraph(),a_vel[a_idir],0,*m_RobinBCData,a_dx[a_idir],a_dit);
    }
  }
  else
  {
    getFaceVel(a_faceFlux,a_face,a_vel,a_probLo,a_dx,a_idir,a_icomp,a_time,a_side,a_doDivFreeOutflow);
  }
}

/*********/

//////called by MAC projection solver for domain bc on phi
void InflowOutflowPoissonDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                               const BaseFab<Real>&  a_phi,
                                               const RealVect&       a_probLo,
                                               const RealVect&       a_dx,
                                               const int&            a_idir,
                                               const Side::LoHiSide& a_side,
                                               const DataIndex&      a_dit,
                                               const Real&           a_time,
                                               const bool&           a_useHomogeneous)
{
  //for phi: outflow  dirichlet. inflow neumann
  bool isOutflow= ((a_side==Side::Hi) && (a_idir==m_flowDir));

  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      // diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
}

//////called by MAC projection solver for domain bc on phi
bool InflowOutflowPoissonDomainBC::
isDirichletDom(const VolIndex&   a_ivof,
               const VolIndex&   a_jvof,
               const EBCellFAB&  a_phi) const
{
  // we don't really need the FaceIndex object, just trying to save some codes here
  FaceIndex face(a_ivof,a_jvof);
  int side = face.faceSign(a_ivof);
  int idir = face.direction();

  bool isOutflow= ((side==1) && (idir==m_flowDir));
  // this logic is valid only for the projection
  if (!isOutflow)
    {
      return false;
    }
  else
    {
      return true;
    }
}

//called by projection solve for irreg x domain
void InflowOutflowPoissonDomainBC::getFaceFlux(Real&                 a_faceFlux,
                                               const VolIndex&       a_vof,
                                               const int&            a_comp,
                                               const EBCellFAB&      a_phi,
                                               const RealVect&       a_probLo,
                                               const RealVect&       a_dx,
                                               const int&            a_idir,
                                               const Side::LoHiSide& a_side,
                                               const DataIndex&      a_dit,
                                               const Real&           a_time,
                                               const bool&           a_useHomogeneous)
{
  bool isOutflow= ((a_side==Side::Hi) && (a_idir==m_flowDir));

  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      // diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
      //                               a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
}

//called by projection solve for irreg x domain
void InflowOutflowPoissonDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
                                                    const VolIndex&       a_vof,
                                                    const int&            a_comp,
                                                    const EBCellFAB&      a_phi,
                                                    const RealVect&       a_probLo,
                                                    const RealVect&       a_dx,
                                                    const int&            a_idir,
                                                    const Side::LoHiSide& a_side,
                                                    const DataIndex&      a_dit,
                                                    const Real&           a_time)
{
  bool isOutflow= ((a_side==Side::Hi) && (a_idir==m_flowDir));

  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit, a_time);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                              a_dx, a_idir, a_side, a_dit, a_time);
      // diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
      //                               a_dx, a_idir, a_side, a_dit, a_time);
    }
}

//called by macEnforceGradientBC, useAreaFrac = false
void InflowOutflowPoissonDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                  const FaceIndex&      a_face,
                                                  const int&            a_comp,
                                                  const EBCellFAB&      a_phi,
                                                  const RealVect&       a_probLo,
                                                  const RealVect&       a_dx,
                                                  const int&            a_idir,
                                                  const Side::LoHiSide& a_side,
                                                  const DataIndex&      a_dit,
                                                  const Real&           a_time,
                                                  const bool&           a_useAreaFrac,
                                                  const RealVect&       a_centroid,
                                                  const bool&           a_useHomogeneous)
{
  bool isOutflow= ((a_side==Side::Hi) && (a_idir==m_flowDir));

  if (!isOutflow)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFaceGradPhi(a_faceFlux, a_face, a_comp,a_phi, a_probLo,
                               a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
  else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.getFaceGradPhi(a_faceFlux, a_face, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useAreaFrac,a_centroid,a_useHomogeneous);
    }
}

//velocity stencil
void InflowOutflowHelmholtzDomainBC::getFluxStencil(      VoFStencil&      a_stencil,
                                                    const VolIndex&        a_vof,
                                                    const int&             a_comp,
                                                    const RealVect&        a_dx,
                                                    const int&             a_idir,
                                                    const Side::LoHiSide&  a_side,
                                                    const EBISBox&         a_ebisBox)
{
  bool isOutflow= ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isSlipWall= ((a_idir!=m_flowDir) && (a_idir != a_comp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));

  if (!m_adjointSolver)
  {
    if ((isOutflow || isSlipWall))
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0);
      neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
    else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.setEBOrder(2);
      diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
  }
  else // adjoint solver
  {
    if (isOutflow || isSlipWall)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.0); // this value doesn't actually matter for Neumann
      neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
    else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      diriBC.setEBOrder(2);
      diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_dx, a_idir, a_side, a_ebisBox);
    }
  }
}

///
/**
   This never gets called.  InflowOutflowPoissonDomainBC::getFaceVel takes care of it.
*/
void
InflowOutflowHelmholtzDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side,
           const bool&           a_doDivFreeOutflow)
{
  MayDay::Error("InflowOutflowHelmholtzDomainBC::getFaceVel is not needed for Helmholtz.");
}
///
/**
   This is called by EBAMRPoissonOp::applyDomainFlux in EBAMRPoissonOp::applyOp
   for reg cells.
   For boundary conditions on velocity in viscous operator.
*/
void InflowOutflowHelmholtzDomainBC::getFaceFlux(BaseFab<Real>&        a_faceFlux,
                                                 const BaseFab<Real>&  a_phi,
                                                 const RealVect&       a_probLo,
                                                 const RealVect&       a_dx,
                                                 const int&            a_idir,
                                                 const Side::LoHiSide& a_side,
                                                 const DataIndex&      a_dit,
                                                 const Real&           a_time,
                                                 const bool&           a_useHomogeneous)
{
  int velcomp =  DirichletPoissonEBBC::s_velComp;

  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((a_idir!=m_flowDir) && (a_idir != velcomp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));

  bool isVelNeum =  ((isOutflow || isSlipWall));

  if (!m_adjointSolver)
  {
    //vel: outflow is neumann. all others dirichlet.  inflow uses inflow vel as the value

/*
  // flip inflow with outflow for adjoint
  if(m_adjointSolver)
  {
    bool tmp = isInflow;
    isInflow = isOutflow;
    isOutflow = tmp;
  }
*/
    if (isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
    }
    else if (isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if (velcomp==m_flowDir)
      {
        if (!m_doPoiseInflow && !m_doWomersleyInflow)
        {
          diriBC.setValue(m_inflowVel);
        }
        else if (m_doPoiseInflow)
        {
          diriBC.setFunction(m_poiseInflowFunc);
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

          diriBC.setValue(m_inflowVel*VelMult/2);
        }
      }
      else
      {
        diriBC.setValue(0.0);
      }

      //basefab flux--EBAMRPoissonOp::applyDomainFlux calls this directly for viscous operator
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
      else
      {
        diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
    }
    else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
      else
      {
        diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
    }
  } // end !m_adjointSolver
  else // m_adjointSolver
  {
    if (isVelNeum)
    {
      // enforcing Robin BC
      if (s_higherOrderHelmBC)
      {
        RobinPoissonBCFunctions::getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous, *m_RobinBCData);
      }
      else
      {
        RobinPoissonBCFunctions::getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous, *m_RobinBCData);
      }
    }
    else
    {  
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
      else
      {
        diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
    }
  }
}

//called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
void InflowOutflowHelmholtzDomainBC::getFaceFlux(Real&                 a_faceFlux,
                                                 const VolIndex&       a_vof,
                                                 const int&            a_comp,
                                                 const EBCellFAB&      a_phi,
                                                 const RealVect&       a_probLo,
                                                 const RealVect&       a_dx,
                                                 const int&            a_idir,
                                                 const Side::LoHiSide& a_side,
                                                 const DataIndex&      a_dit,
                                                 const Real&           a_time,
                                                 const bool&           a_useHomogeneous)
{
  int velcomp =  DirichletPoissonEBBC::s_velComp;

  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((a_idir!=m_flowDir) && (a_idir != velcomp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));

  bool isVelNeum =  ((isOutflow || isSlipWall));

  if (!m_adjointSolver)
  {
    //vel: outflow is Neumann. all others Dirichlet

/*
  // flip inflow with outflow for adjoint
  if(m_adjointSolver)
  {
    bool tmp = isInflow;
    isInflow = isOutflow;
    isOutflow = tmp;
  }
*/

    if (isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    }
    else if (isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if (velcomp==m_flowDir)
      {
        if (!m_doPoiseInflow && !m_doWomersleyInflow)
        {
          diriBC.setValue(m_inflowVel);
        }
        else if (m_doPoiseInflow)
        {
          diriBC.setFunction(m_poiseInflowFunc);
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

          diriBC.setValue(m_inflowVel*VelMult/2);
        }
      }
      else
      {
        diriBC.setValue(0.0);
      }
      //called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
      else
      {
        diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
      }
    }
    else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
      else
      {
        diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
      }
    }
  } // end !adjointSolver
  else
  {
    if (isVelNeum)
    {
      if (s_higherOrderHelmBC)
      {
        RobinPoissonBCFunctions::getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous, *m_RobinBCData);
      }
      else
      {
        RobinPoissonBCFunctions::getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous, *m_RobinBCData);
      }
    }
    else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, a_useHomogeneous);
      }
      else
      {
        diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
      }
    }
  }
}

//////called by Helm solver for domain bc on phi
bool InflowOutflowHelmholtzDomainBC::
isDirichletDom(const VolIndex&   a_ivof,
               const VolIndex&   a_jvof,
               const EBCellFAB&  a_phi) const
{

  FaceIndex face(a_ivof,a_jvof);
  int side = face.faceSign(a_ivof);
  int idir = face.direction();
  int velcomp =  DirichletPoissonEBBC::s_velComp;
  bool isOutflow = ((side== 1) && (idir==m_flowDir));
  bool isSlipWall = ((idir!=m_flowDir) && (idir != velcomp) && ((m_doSlipWallsHi[idir]==1 && side == 1)||(m_doSlipWallsLo[idir]==1 && side == -1)));

  bool isVelNeum =  ((isOutflow || isSlipWall));

  if (!m_adjointSolver)
  {

    if (isVelNeum)
    {
      return false;
    }
    else
    {
      return true;
    }
  }
  else
  {
    if (isVelNeum)
    {
      return true;
    }
    else
    {
      return true;
    }
  }
}

//called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
void InflowOutflowHelmholtzDomainBC::getInhomFaceFlux(Real&                 a_faceFlux,
                                                      const VolIndex&       a_vof,
                                                      const int&            a_comp,
                                                      const EBCellFAB&      a_phi,
                                                      const RealVect&       a_probLo,
                                                      const RealVect&       a_dx,
                                                      const int&            a_idir,
                                                      const Side::LoHiSide& a_side,
                                                      const DataIndex&      a_dit,
                                                      const Real&           a_time)
{
  int velcomp =  DirichletPoissonEBBC::s_velComp;

  bool isOutflow =  ((a_side==Side::Hi) && (a_idir==m_flowDir));
  bool isInflow =  ((a_side==Side::Lo) && (a_idir==m_flowDir));
  bool isSlipWall = ((a_idir!=m_flowDir) && (a_idir != velcomp) && ((m_doSlipWallsHi[a_idir]==1 && a_side == Side::Hi)||(m_doSlipWallsLo[a_idir]==1 && a_side == Side::Lo)));

  bool isVelNeum =  ((isOutflow || isSlipWall));

  if (!m_adjointSolver)
  {
    //vel: outflow is Neumann. all others Dirichlet

/*
  // flip in case of adjoint
  if (m_adjointSolver)
  {
    bool tmp = isInflow;
    isInflow = isOutflow;
    isOutflow = tmp;
  }
*/

    if (isVelNeum)
    {
      NeumannPoissonDomainBC neumannBC;
      neumannBC.setValue(0.);
      neumannBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                                 a_dx, a_idir, a_side, a_dit,a_time);
    }
    else if (isInflow)
    {
      DirichletPoissonDomainBC diriBC;
      if (velcomp==m_flowDir)
      {
        if (!m_doPoiseInflow && !m_doWomersleyInflow)
        {
          diriBC.setValue(m_inflowVel);
        }
        else if (m_doPoiseInflow)
        {
          diriBC.setFunction(m_poiseInflowFunc);
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

          diriBC.setValue(m_inflowVel*VelMult/2);
        }
      }
      else
      {
        diriBC.setValue(0.0);
      }
      //called by EBAMRPoissonOp::applyOp for viscous Helmholtz EB x domain
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time);
      }
      else
      {
        diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time);
      }
    }
    else
    {
      //wall bc no slip
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time);
      }
      else
      {
        diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time);
      }
    }
  }  // end !adjointSolver
  else
  {
    if (isVelNeum)
    {
      if (s_higherOrderHelmBC)
      {
        RobinPoissonBCFunctions::getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, *m_RobinBCData);
      }
      else
      {
        RobinPoissonBCFunctions::getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time, *m_RobinBCData);
      }
    }
    else
    {
      DirichletPoissonDomainBC diriBC;
      diriBC.setValue(0.0);
      if (s_higherOrderHelmBC)
      {
        diriBC.getHigherOrderInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit, a_time);
      }
      else
      {
        diriBC.getInhomFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time);
      }
    }
  }
}

void InflowOutflowHelmholtzDomainBC::getFaceGradPhi(Real&                 a_faceFlux,
                                                    const FaceIndex&      a_face,
                                                    const int&            a_comp,
                                                    const EBCellFAB&      a_phi,
                                                    const RealVect&       a_probLo,
                                                    const RealVect&       a_dx,
                                                    const int&            a_idir,
                                                    const Side::LoHiSide& a_side,
                                                    const DataIndex&      a_dit,
                                                    const Real&           a_time,
                                                    const bool&           a_useAreaFrac,
                                                    const RealVect&       a_centroid,
                                                    const bool&           a_useHomogeneous)
{
  MayDay::Error("InflowOutflowHelmholtzDomainBC::getFaceGradPhi not needed");
}

#include "NamespaceFooter.H"
