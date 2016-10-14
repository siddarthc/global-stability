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

#ifndef  _InflowOutflowIBC_H__
#define  _InflowOutflowIBC_H__

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"
#include "EBPhysIBC.H"
#include "InflowOutflowParams.H"

#include "NamespaceHeader.H"

class InflowOutflowIBC : public EBPhysIBC
{
public:
  virtual ~InflowOutflowIBC();

  InflowOutflowIBC(const InflowOutflowParams& a_params);

    ///
  void define(const ProblemDomain&  a_domain,
              const RealVect&       a_dx);

  ///  For every box in this level, this function is called.
  void fluxBC(EBFluxFAB&            a_flux,
              const EBCellFAB&      a_Wcenter,
              const EBCellFAB&      a_Wextrap,
              const Side::LoHiSide& a_sd,
              const Real&           a_time,
              const EBISBox&        a_ebisBox,
              const DataIndex&      a_dit,
              const Box&            a_box,
              const Box&            a_faceBox,
              const int&            a_dir);

  /// Initialize
  void initialize(LevelData<EBCellFAB>& a_conState,
                  const EBISLayout& a_ebisl) const;

  /// Set boundary slopes
  /**
     The boundary slopes in a_dW are already set to one sided difference
     approximations.  If this function doesn't change them they will be
     used for the slopes at the boundaries.
  */
  void setBndrySlopes(EBCellFAB&       a_deltaPrim,
                      const EBCellFAB& a_primState,
                      const EBISBox&   a_ebisBox,
                      const Box&       a_box,
                      const int&       a_dir);

protected:
  InflowOutflowParams  m_params;

  bool m_isDefined;

  Real m_dx;
  ProblemDomain  m_domain;

private:
  InflowOutflowIBC()
  {
    MayDay::Error("invalid operator");
  }

  //disallowed for all the usual reasons
  void operator=(const InflowOutflowIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }
  InflowOutflowIBC(const InflowOutflowIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#include "NamespaceFooter.H"
#endif
