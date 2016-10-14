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

#include "InflowOutflowVelocityDomBC.H"
#include "BoxIterator.H"
#include "NeumannViscousTensorDomainBC.H"
#include "DirichletViscousTensorDomainBC.H"
#include "NeumannViscousTensorDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletViscousTensorEBBC.H"
#include "NeumannViscousTensorEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"
#include "SlipWallViscousTensorDomainBC.H"

#include "NamespaceHeader.H"

class InflowFunction : public BaseBCFuncEval
{
public:
  InflowFunction(const RealVect& a_inflowVel)
  {
    m_inflowVel = a_inflowVel;
  }

  virtual Real value(const RealVect& a_point,
                     const int&      a_comp) const
  {
    Real retval = m_inflowVel[a_comp];
//    pout() << "set inflow vel = " << retval << " for comp = " << a_comp << endl;
    return retval;
  }

  virtual Real derivative(const RealVect& a_point,
                          const int&      a_comp,
                          const int&      a_derivDir) const
  {
    return 0.;
  }

protected:

  RealVect m_inflowVel;
};

void InflowOutflowVelocityDomBC::
whereAMI(bool& a_atInflow,
         bool& a_atOutflow,
         bool& a_atWall,
         bool& a_atSlipWall,
         const int& a_dir,
         const Side::LoHiSide& a_side)
{
  //figure out if we are on a wall and if its an inflow and the inflow value
  a_atInflow     = ((a_side == Side::Lo) && (a_dir==m_params.m_inflowDir));
  a_atOutflow    = ((a_side == Side::Hi) && (a_dir==m_params.m_outflowDir));
  a_atWall       = ((!(a_atInflow)) && (!(a_atOutflow)) && (((a_side == Side::Lo) && (m_params.m_wallBCLo[a_dir]==1)) || ((a_side == Side::Hi) && (m_params.m_wallBCHi[a_dir]==1))));
  a_atSlipWall   = ((a_atWall) && (((a_side == Side::Lo) && (m_params.m_slipWallLo[a_dir]==1)) || ((a_side == Side::Hi) && (m_params.m_slipWallHi[a_dir]==1))));

}

void InflowOutflowVelocityDomBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  bool atInflow, atOutflow, atWall, atSlipWall;
  whereAMI(atInflow, atOutflow, atWall, atSlipWall, a_idir, a_side); 
  if (atInflow)
    {
//      pout() << "at inflow boundary for dir = " << a_idir << " side = "<< a_side << endl;
      DirichletViscousTensorDomainBC diriBC;
      RefCountedPtr<BaseBCFuncEval> funk(new InflowFunction (m_params.m_inflowVel));
      diriBC.setFunction(funk);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    } // end atInflow
  else if (atOutflow)
    {
//      pout() << "at outflow boundary for dir = " << a_idir << "side = " << a_side << endl; 
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    } // end atOutflow
  else if (atWall)
    {
      if (atSlipWall)
        {
//          pout() << "at slip bndry for dir = " << a_idir << "side = " << a_side << endl;
          NeumannViscousTensorDomainBC  neumBC;
          Real value = 0;
          neumBC.setValue(value);

          neumBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          neumBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        } // end atSlipWall
      else // noSlip
        {
//          pout() << "ar no slip bndry for dir = " << a_idir << "side = " << a_side << endl;
          DirichletViscousTensorDomainBC diriBC;
          Real value = 0;
          diriBC.setValue(value);

          diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        } // end atNoSlipWall
    } // end atWall
  else
    {
      pout() << "invalid BC in InflowOutflowVelocityDomainBC" << endl;
    }
}

void InflowOutflowVelocityDomBC::
getFaceFlux(Real&                 a_faceFlux,
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
  bool atInflow, atOutflow, atWall, atSlipWall;
  whereAMI(atInflow, atOutflow, atWall, atSlipWall, a_idir, a_side);
  if (atInflow)
    {
      DirichletViscousTensorDomainBC diriBC;
      RefCountedPtr<BaseBCFuncEval> funk(new InflowFunction (m_params.m_inflowVel));
      diriBC.setFunction(funk);

      diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                         a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    } // end atInflow
  else if (atOutflow)
    { 
      NeumannViscousTensorDomainBC neumannBC;
      neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
      neumannBC.setValue(0.0);
      neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                            a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
    } // end atOutflow
  else if (atWall)
    {
      if (atSlipWall)
        {
/*
          SlipWallViscousTensorDomainBC slipbc;
          slipbc.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          slipbc.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                             a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
*/
          NeumannViscousTensorDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          neumannBC.setValue(0.0);
          neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                                a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        } // end atSlipWall
      else // noSlip
        {
          DirichletViscousTensorDomainBC diriBC;
          Real value = 0;
          diriBC.setValue(value);

          diriBC.setCoef(m_eblg, m_beta, m_eta, m_lambda);
          diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                             a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
        } // end atNoSlipWall
    } // end atWall
  else
    {
      pout() << "invalid BC in InflowOutflowVelocityDomainBC" << endl;
    }
}

#include "NamespaceFooter.H"
