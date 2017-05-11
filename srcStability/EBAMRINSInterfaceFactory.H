/*
 *
 *
 *
 *
 */

#ifndef _EBAMRINS_Interface_Factory_H_
#define _EBAMRINS_Interface_Factory_H_

#include "ChomboSolverInterfaceFactory.H"
#include "EBAMRINSInterface.H"

class EBAMRINSInterfaceFactory : public ChomboSolverInterfaceFactory
{
public:

  EBAMRINSInterfaceFactory(const AMRParameters& a_params,
                           const RefCountedPtr<EBIBCFactory>& a_baseflowIBC,
                           const RefCountedPtr<EBIBCFactory>& a_solverIBC,
                           const ProblemDomain& a_coarsestDomain,
                           Real                 a_viscosity,
                           int                  a_snapshotInterval,
                           bool                 a_doLinINS,
                           bool                 a_doAdjoint,
                           bool                 a_doTransientGrowth,
                           bool                 a_doFirstOrderFreDeriv = false,
                           const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance())

  {
    m_params          = a_params;
    m_baseflowIBCFact = a_baseflowIBC;
    m_solverIBCFact   = a_solverIBC;
    m_coarsestDomain  = a_coarsestDomain;
    m_viscosity       = a_viscosity;
    m_ebisPtr         = a_ebisPtr;
    m_doLinINS        = a_doLinINS;
    m_doAdjoint       = a_doAdjoint;
    m_doTransientGrowth = a_doTransientGrowth;
    m_doFirstOrderFreDeriv = a_doFirstOrderFreDeriv;

    m_snapshotInterval = a_snapshotInterval;
  }

/*********/

  ~EBAMRINSInterfaceFactory()
  {;}

/*********/

  ChomboSolverInterface* create() const
  {
    EBAMRINSInterface* retval = new EBAMRINSInterface(m_params,
                                                      m_baseflowIBCFact,
                                                      m_solverIBCFact,
                                                      m_coarsestDomain,
                                                      m_viscosity,
                                                      m_snapshotInterval,
                                                      m_doLinINS,
                                                      m_doAdjoint,
                                                      m_doTransientGrowth,
                                                      m_doFirstOrderFreDeriv,
                                                      m_ebisPtr);

    return static_cast<ChomboSolverInterface*>(retval);
  }

private:

  AMRParameters               m_params;
  RefCountedPtr<EBIBCFactory> m_baseflowIBCFact;
  RefCountedPtr<EBIBCFactory> m_solverIBCFact;
  ProblemDomain               m_coarsestDomain;
  Real                        m_viscosity;
  int                         m_snapshotInterval;
  bool                        m_doLinINS;
  bool                        m_doAdjoint;
  bool                        m_doTransientGrowth;
  bool                        m_doFirstOrderFreDeriv;
  const EBIndexSpace*         m_ebisPtr;

  EBAMRINSInterfaceFactory()
  {
    MayDay::Error("EBAMRINSInterfaceFactory::Invalid Constructor");
  }
};

#endif
