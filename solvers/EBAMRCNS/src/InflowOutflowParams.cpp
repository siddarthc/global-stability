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
  pp.get("gamma"           ,   a_params.m_gamma);
  pp.get("inflowDir"       ,   a_params.m_inflowDir);
  pp.get("outflowDir"      ,   a_params.m_outflowDir);
  pp.get("inflowDense"     ,   a_params.m_inflowDense);
  pp.get("inflowPress"     ,   a_params.m_inflowPress);

  a_params.m_outflowPress = a_params.m_inflowPress;
  pp.query("outflowPress" , a_params.m_outflowPress);

  a_params.m_totPresBC = false;
  pp.query("totalPresBC"     ,    a_params.m_totPresBC);

  vector<Real> inflowVelocity(SpaceDim);
  vector<Real> initVelocity(SpaceDim);

  pp.getarr("inflowVel", inflowVelocity, 0, SpaceDim);
  pp.getarr("initVel", initVelocity, 0, SpaceDim);
  for(int idir = 0; idir < SpaceDim; idir++)
  {
      a_params.m_inflowVel[idir] = inflowVelocity[idir];
      a_params.m_initVel[idir] = initVelocity[idir];
  }

  pp.get("initPress"       ,   a_params.m_initPress);
  pp.get("initDense"       ,   a_params.m_initDense);
  pp.get("specific_heat"    ,   a_params.m_specificHeat);

  Real R = a_params.m_specificHeat*(a_params.m_gamma - 1.0);
  Real denom0 = a_params.m_initDense*R;
  Real denom1 = a_params.m_inflowDense*R;
  a_params.m_initTemp = a_params.m_initPress/denom0;
  a_params.m_inflowTemp = a_params.m_inflowPress/denom1;

  vector<int> wallLow(SpaceDim);
  vector<int> wallHigh(SpaceDim);
  vector<int> slipWallLow(SpaceDim);
  vector<int> slipWallHigh(SpaceDim);
  vector<int> adiabaticBCLow(SpaceDim);
  vector<int> adiabaticBCHigh(SpaceDim);
  vector<int> isothermalBCLow(SpaceDim);
  vector<int> isothermalBCHigh(SpaceDim);
  vector<Real> isoTempBCLow(SpaceDim);
  vector<Real> isoTempBCHigh(SpaceDim);
  vector<Real> domLength(SpaceDim);

  pp.getarr("wallBCLo", wallLow, 0, SpaceDim);
  pp.getarr("wallBCHi", wallHigh, 0, SpaceDim);
  pp.getarr("slipWallLo", slipWallLow, 0, SpaceDim);
  pp.getarr("slipWallHi", slipWallHigh, 0, SpaceDim);
  pp.getarr("do_adiabaticBCLo", adiabaticBCLow, 0, SpaceDim);
  pp.getarr("do_adiabaticBCHi", adiabaticBCHigh, 0, SpaceDim);
  pp.getarr("do_isothermalBCLo", isothermalBCLow, 0, SpaceDim);
  pp.getarr("do_isothermalBCHi", isothermalBCHigh, 0, SpaceDim);
  pp.getarr("isoTempBCLo", isoTempBCLow, 0, SpaceDim);
  pp.getarr("isoTempBCHi", isoTempBCHigh, 0, SpaceDim);
  pp.getarr("domain_length", domLength, 0, SpaceDim);

  for (int idir = 0; idir < SpaceDim; idir++)
  {
      a_params.m_wallBCLo[idir] = wallLow[idir];
      a_params.m_wallBCHi[idir] = wallHigh[idir];
      a_params.m_slipWallLo[idir] = slipWallLow[idir];
      a_params.m_slipWallHi[idir] = slipWallHigh[idir];
      a_params.m_adiabaticBCLo[idir] = adiabaticBCLow[idir];
      a_params.m_adiabaticBCHi[idir] = adiabaticBCHigh[idir];
      a_params.m_isothermalBCLo[idir] = isothermalBCLow[idir];
      a_params.m_isothermalBCHi[idir] = isothermalBCHigh[idir];
      a_params.m_isoTempBCLo[idir] = isoTempBCLow[idir];
      a_params.m_isoTempBCHi[idir] = isoTempBCHigh[idir];
      a_params.m_domLength[idir] = domLength[idir];
  }
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
