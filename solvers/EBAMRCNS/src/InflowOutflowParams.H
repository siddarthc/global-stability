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
    m_gamma = 1.4;
    m_inflowDir = 0;
    m_outflowDir =0;
    m_inflowDense = 1.225;
    m_inflowPress = 101325;
    m_outflowPress = 101325;
    m_inflowTemp = 300;
    for (int i = 0; i < SpaceDim; i++)
      {
        m_inflowVel[i] = 0;
        m_initVel[i] = 0;
        m_slipWallLo[i] = 0;
        m_slipWallHi[i] = 0;
        m_adiabaticBCLo[i] = 0;
        m_adiabaticBCHi[i] = 0;
        m_isothermalBCLo[i] = 0;
        m_isothermalBCHi[i] = 0;
        m_isoTempBCLo[i] = 0;
        m_isoTempBCHi[i] = 0;
        if (i <= 1)
          {
            m_wallBCHi[i] = i;
            m_wallBCLo[i] = i;
          }
        else
          {
            m_wallBCHi[i] = i-1;
            m_wallBCLo[i] = i-1;
          }
      }
    m_initTemp = 300;
    m_initPress = 101325;
    m_initDense = 1.225;
    m_specificHeat = 717.5;
    m_totPresBC = 0;
    m_verbosity = 1;
    m_doInitWhiteNoise = 0;
    m_initSigma = 0;
    m_doInflowWhiteNoise = 0;
    m_inflowSigma = 0;
  }

  ~InflowOutflowParams(){}

  Real m_gamma;
  int  m_inflowDir;
  int  m_outflowDir;
  Real m_inflowDense;
  Real m_inflowPress;
  Real m_outflowPress;
  Real m_inflowTemp;
  RealVect m_inflowVel;
  RealVect m_initVel;
  Real m_initTemp;
  Real m_initPress;
  Real m_initDense;
  Real m_specificHeat;
  IntVect m_wallBCHi;
  IntVect m_wallBCLo;
  IntVect m_slipWallLo;
  IntVect m_slipWallHi;
  IntVect m_adiabaticBCLo;
  IntVect m_adiabaticBCHi;
  IntVect m_isothermalBCLo;
  IntVect m_isothermalBCHi;
  RealVect m_isoTempBCLo;
  RealVect m_isoTempBCHi;
  RealVect m_domLength;
  bool m_totPresBC;
  int m_verbosity;
  bool m_doInflowWhiteNoise;
  Real m_inflowSigma;
  bool m_doInitWhiteNoise;
  Real m_initSigma;
};

///get stuff from input file
extern void ParseInflowOutflowParams(InflowOutflowParams& a_params);

#include "NamespaceFooter.H"
#endif
