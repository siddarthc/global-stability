/*
 * This pig is added to protect this code from the bugs.
 *      To make best use of this pig, compile with opt and fPIG flags
 *                               _
 *       _._ _..._ .-',     _.._(`))
 *      '-. `     '  /-._.-'    ',/
 *         )         \            '.
 *        / _    _    |             \
 *       |  o    o    /              |
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

#ifndef _InflowOutflowTemperatureDomBC_H_
#define _InflowOutflowTemperatureDomBC_H_

#include "BaseDomainBC.H"
#include "RefCountedPtr.H"
#include "InflowOutflowParams.H"

#include "NamespaceHeader.H"

class InflowOutflowTemperatureDomBC : public ConductivityBaseDomainBC
{
public:

  virtual ~InflowOutflowTemperatureDomBC()
  {;}

  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

///
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
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
    //this is for projections and such.
    //not used here.
    MayDay::Error();
  }

  virtual void getFaceVel(Real&                 a_faceFlux,
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
    //this is for projections and such.
    //not used here.
    MayDay::Error();
  }

  InflowOutflowTemperatureDomBC(const InflowOutflowParams& a_params)
  {
    m_params = a_params;
  }

private:

  void
  whereAMI(bool& a_atInflow,
           bool& a_atOutflow,
           bool& a_atAdiabaticBndry,
           bool& a_atIsothermalBndry,
           const int&            a_idir,
           const Side::LoHiSide& a_side);

  InflowOutflowTemperatureDomBC(){;};
  InflowOutflowParams m_params;
};

class InflowOutflowTemperatureDomBCFactory : public BaseDomainBCFactory
{
public:

  InflowOutflowTemperatureDomBCFactory(const InflowOutflowParams& a_params)
  {
    m_params = a_params;
  }

  virtual ~InflowOutflowTemperatureDomBCFactory()
  {
  }

  virtual InflowOutflowTemperatureDomBC* create(const ProblemDomain& a_domain,
                                                const EBISLayout&    a_layout,
                                                const RealVect&      a_dx)
  {
    InflowOutflowTemperatureDomBC* retval = new InflowOutflowTemperatureDomBC(m_params);
    return retval;
  }

private:

  InflowOutflowTemperatureDomBCFactory(){;}
  InflowOutflowParams m_params;
};

#include "NamespaceFooter.H"
#endif
