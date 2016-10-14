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

#include "InflowOutflowTemperatureDomBC.H"
#include "BoxIterator.H"
#include "NeumannConductivityDomainBC.H"
#include "DirichletConductivityDomainBC.H"
#include "EBArith.H"
#include "Stencils.H"
#include "DirichletConductivityEBBC.H"
#include "VoFIterator.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

void InflowOutflowTemperatureDomBC::
whereAMI(bool& a_atInflow,
         bool& a_atOutflow,
         bool& a_atAdiabaticBndry,
         bool& a_atIsothermalBndry,
         const int& a_dir,
         const Side::LoHiSide& a_side)
{
  //figure out if we are on a wall and if its an inflow and the inflow value
  a_atInflow     = ((a_side == Side::Lo) && (a_dir==m_params.m_inflowDir));
  a_atOutflow    = ((a_side == Side::Hi) && (a_dir==m_params.m_outflowDir));
  a_atAdiabaticBndry = ((!(a_atInflow)) && (!(a_atOutflow)) && (((a_side == Side::Lo) && (m_params.m_adiabaticBCLo[a_dir]==1)) || ((a_side == Side::Hi) && (m_params.m_adiabaticBCHi[a_dir]==1))));
  a_atIsothermalBndry = ((!(a_atInflow)) && (!(a_atOutflow)) && (!(a_atAdiabaticBndry)) && (((a_side == Side::Lo) && (m_params.m_isothermalBCLo[a_dir]==1)) || ((a_side == Side::Hi) && (m_params.m_isothermalBCHi[a_dir]==1))));

}

void InflowOutflowTemperatureDomBC::
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
  bool atInflow, atOutflow, atAdiabaticBndry, atIsothermalBndry;
  whereAMI(atInflow, atOutflow, atAdiabaticBndry, atIsothermalBndry, a_idir, a_side);
  if (atInflow)
  {
/*
    DirichletConductivityDomainBC diriBC;
    Real value = m_params.m_inflowTemp;
    diriBC.setValue(value);
    diriBC.setCoef(m_eblg, m_beta, m_bcoef);
    diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous); 
*/
    NeumannConductivityDomainBC neumannBC;
    neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
    neumannBC.setValue(0.0);
    neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  } // end atInflow
  else if (atOutflow)
  {
    NeumannConductivityDomainBC neumannBC;
    neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
    neumannBC.setValue(0.0);
    neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  }  // end atOutflow
  else if (atAdiabaticBndry)
  {
    NeumannConductivityDomainBC neumannBC;
    neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
    neumannBC.setValue(0.0);
    neumannBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  } // end atAdiabaticBndry
  else if (atIsothermalBndry)
  {
    DirichletConductivityDomainBC diriBC;
    Real value = (a_side == Side::Lo) ? m_params.m_isoTempBCLo[a_idir] : m_params.m_isoTempBCHi[a_idir];
    diriBC.setValue(value);
    diriBC.setCoef(m_eblg, m_beta, m_bcoef);
    diriBC.getFaceFlux(a_faceFlux, a_phi, a_probLo, a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  } // end atIsothermalBndry
  else 
  {
    pout() << "Invalid BC in InflowOutflowTemperatureDomBC" << endl;
  }
}

void InflowOutflowTemperatureDomBC::
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
  bool atInflow, atOutflow, atAdiabaticBndry, atIsothermalBndry;
  whereAMI(atInflow, atOutflow, atAdiabaticBndry, atIsothermalBndry, a_idir, a_side);
  if (atInflow)
  {

    DirichletConductivityDomainBC diriBC;
    Real value = m_params.m_inflowTemp;
    diriBC.setValue(value);
    diriBC.setCoef(m_eblg, m_beta, m_bcoef);
    diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                       a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
/*
    NeumannConductivityDomainBC neumannBC;
    neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
    neumannBC.setValue(0.0);
    neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                          a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
*/
  } // end atInflow
  else if (atOutflow)
  {
    NeumannConductivityDomainBC neumannBC;
    neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
    neumannBC.setValue(0.0);
    neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                          a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  }  // end atOutflow
  else if (atAdiabaticBndry)
  {
    NeumannConductivityDomainBC neumannBC;
    neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
    neumannBC.setValue(0.0);
    neumannBC.getFaceFlux(a_faceFlux, a_vof, a_comp,a_phi, a_probLo,
                          a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  } // end atAdiabaticBndry
  else if (atIsothermalBndry)
  {
    DirichletConductivityDomainBC diriBC;
    Real value = (a_side == Side::Lo) ? m_params.m_isoTempBCLo[a_idir] : m_params.m_isoTempBCHi[a_idir];
    diriBC.setValue(value);
    diriBC.setCoef(m_eblg, m_beta, m_bcoef);
    diriBC.getFaceFlux(a_faceFlux, a_vof, a_comp, a_phi, a_probLo,
                       a_dx, a_idir, a_side, a_dit,a_time,a_useHomogeneous);
  } // end atIsothermalBndry
  else
  {
    pout() << "Invalid BC in InflowOutflowTemperatureDomBC" << endl;
  }
}

#include "NamespaceFooter.H"
