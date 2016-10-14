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

#ifndef _InflowOutflowParams_H_
#define _InflowOutflowParams_H_

#include "Vector.H"
#include "RealVect.H"
#include "REAL.H"
#include "IntVect.H"
#include "NamespaceHeader.H"

class InflowOutflowParams
{
public:

  InflowOutflowParams()
  {
    // CAUTION :- Always parse inflow-outflow params to avoid
    // an error that might result from using default parameters !
    m_doInitWhiteNoise = 0;
    m_initSigma = 0;
    m_doInflowWhiteNoise = 0;
    m_inflowSigma = 0;
  }

  bool m_doInflowWhiteNoise;
  Real m_inflowSigma;
  bool m_doInitWhiteNoise;
  Real m_initSigma;

};

///get stuff from input file
extern void ParseInflowOutflowParams(InflowOutflowParams& a_params);

#include "NamespaceFooter.H"
#endif
