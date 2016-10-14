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

#include "InflowOutflowIBCFactory.H"
#include "InflowOutflowIBC.H"

#include "NamespaceHeader.H"

/***************/
InflowOutflowIBCFactory::
InflowOutflowIBCFactory(const InflowOutflowParams& a_params) 
   :EBPhysIBCFactory()
{
  m_params = a_params;
}
/****************/
InflowOutflowIBCFactory::
~InflowOutflowIBCFactory()
{;}
/****************/
EBPhysIBC*
InflowOutflowIBCFactory::
create() const
{
  InflowOutflowIBC* retval = new InflowOutflowIBC(m_params);
  return static_cast<EBPhysIBC*>(retval);
}
/****************/

#include "NamespaceFooter.H"
