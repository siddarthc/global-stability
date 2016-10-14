/* 
 * This pig is added to protect this code from the bugs.
 * 	To make best use of this pig, compile with opt and fPIG flags	 
 *	                         _
 *	 _._ _..._ .-',     _.._(`))
 *	'-. `     '  /-._.-'    ',/
 *	   )         \            '.
 *	  / _    _    |             \
 *	 |  o    o    /              |
 *	 \   .-.                     ;  
 *	  '-('' ).-'       ,'       ;
 *	     '-;           |      .'
 *	        \           \    /
 *	        | 7  .__  _.-\   \
 *	        | |  |  ``/  /`  /
 *	       /,_|  |   /,_/   /
 *	          /,_/      '`-'
 *
 *
 * refer = curious ? 
 * 	   https://gist.github.com/crittelmeyer/b4c8ec6f02c024798ccb 
 * 	   : NULL;
**/

#ifndef _InflowOutflowIBCFACTORY_H
#define _InflowOutflowIBCFACTORY_H

#include "EBCellFAB.H"
#include "EBFaceFAB.H"
#include "REAL.H"
#include "LevelData.H"
#include "InflowOutflowParams.H"
#include "ProblemDomain.H"
#include "EBPhysIBC.H"
#include "InflowOutflowIBC.H"
#include "EBPhysIBCFactory.H"

#include "NamespaceHeader.H"

///
/**
 */
class InflowOutflowIBCFactory: public EBPhysIBCFactory
{
public:

  ///
  ~InflowOutflowIBCFactory();

  ///
  InflowOutflowIBCFactory(const InflowOutflowParams&  a_params);


  ///
  EBPhysIBC* create() const;
protected:
  InflowOutflowParams  m_params;

private:
  InflowOutflowIBCFactory()
  {
    MayDay::Error("Invalid operator");
  }

  void operator=(const InflowOutflowIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }

  InflowOutflowIBCFactory(const InflowOutflowIBCFactory& a_input)
  {
    MayDay::Error("Invalid operator");
  }


};

#include "NamespaceFooter.H"
#endif
