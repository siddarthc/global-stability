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
                           const RefCountedPtr<EBIBCFactory>&  a_IBC,
                           const ProblemDomain& a_coarsestDomain,
                           Real                 a_viscosity,
                           bool                 a_plotSnapshots,
                           bool                 a_doFirstOrderFreDeriv,
                           const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance())

  {
    m_params         = a_params;
    m_ibcFact        = a_IBC;
    m_coarsestDomain = a_coarsestDomain;
    m_viscosity      = a_viscosity;
    m_ebisPtr        = a_ebisPtr;
    m_plotSnapshots  = a_plotSnapshots;
    m_doFirstOrderFreDeriv = a_doFirstOrderFreDeriv;
  }

/*********/

  ~EBAMRINSInterfaceFactory()
  {;}

/*********/

  ChomboSolverInterface* create() const
  {
    EBAMRINSInterface* retval = new EBAMRINSInterface(m_params,
                                                      m_ibcFact,
                                                      m_coarsestDomain,
                                                      m_viscosity,
                                                      m_plotSnapshots,
                                                      m_doFirstOrderFreDeriv,
                                                      m_ebisPtr);

    return static_cast<ChomboSolverInterface*>(retval);
  }

private:

  AMRParameters               m_params;
  RefCountedPtr<EBIBCFactory> m_ibcFact;
  ProblemDomain               m_coarsestDomain;
  Real                        m_viscosity;
  bool                        m_plotSnapshots;
  bool                        m_doFirstOrderFreDeriv;
  const EBIndexSpace*         m_ebisPtr;

  EBAMRINSInterfaceFactory()
  {
    MayDay::Error("EBAMRINSInterfaceFactory::Invalid Constructor");
  }
};

#endif
