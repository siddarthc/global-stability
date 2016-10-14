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

#include "InflowOutflowParams.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

#define PI 3.14159265358979323846  /*   pi  */

void
ParseInflowOutflowParams(InflowOutflowParams& a_params)
{
  ParmParse pp;
  // read params for white noise:-
  pp.get("do_initWhiteNoise", a_params.m_doInitWhiteNoise);
  pp.get("do_inflowWhiteNoise", a_params.m_doInflowWhiteNoise);

  Real initAmp, inflowAmp;
  pp.get("initWhiteNoise", initAmp);
  pp.get("inflowWhiteNoise", inflowAmp);

  a_params.m_initSigma = (SpaceDim == 2) ? initAmp/(sqrt(2.)*3.) : initAmp/(sqrt(3.)*3.);
  a_params.m_inflowSigma = (SpaceDim == 2) ? inflowAmp/(sqrt(2.)*3.) : inflowAmp/(sqrt(3.)*3.);
}

#include "NamespaceFooter.H"
