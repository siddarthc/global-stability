/*
 *
 *
 *
 *
 */

#include "ChomboSFDInterface.H"
#include "EBCellFactory.H"
#include "ChomboSFDOpFunctions.H"
#include "EBLevelDataOps.H"
#include "ParmParse.H"
#include <sstream>

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/

ChomboSFDInterface::ChomboSFDInterface()
{
  m_isDefined = false;
  m_doPIControl = false;
  m_nResetSteps = 0;
}
/*********/
ChomboSFDInterface::~ChomboSFDInterface()
{

  for (int i = 0; i< m_qBar.size(); i++)
    {
      if (m_qBar[i] != NULL) delete m_qBar[i];
      // for PI control
      if (m_qdiffOld[i] != NULL) delete m_qdiffOld[i];
      if (m_qdiffNew[i] != NULL) delete m_qdiffNew[i];
      if (m_qdiffNew[i] != NULL) delete m_qdiffSum[i];
      if (m_qBar_tOld[i] != NULL) delete m_qBar_tOld[i];
      if (m_qdiffSum_tOld[i] != NULL) delete m_qdiffSum_tOld[i];
    }

}
/*********/
void ChomboSFDInterface::
define(int    a_nFilters,
       double a_smallestFilter,
       double a_largestFilter,
       double a_controlCoef,
       int    a_nComp,
       bool   a_doPIControl,
       Vector<Real> a_integralCoef)
{
  m_qBar.resize(a_nFilters);
  m_SFDOp.define(a_nFilters, a_smallestFilter, a_largestFilter, a_controlCoef);
  m_nComp = a_nComp;

  // for PI control
  m_integralCoef = a_integralCoef;
  m_doPIControl = a_doPIControl;
  if (m_doPIControl) 
  {
    m_qdiffOld.resize(a_nFilters);
    m_qdiffNew.resize(a_nFilters);
    m_qdiffSum.resize(a_nFilters);
    m_qBar_tOld.resize(a_nFilters);
    m_qdiffSum_tOld.resize(a_nFilters);
  }
  
  m_isDefined = true;
}
/*********/
void ChomboSFDInterface::
initialize(const LevelData<EBCellFAB>& a_data,
           EBLevelGrid&                a_eblg)
{
  CH_assert(m_isDefined)

  EBCellFactory fact(a_eblg.getEBISL());
  const IntVect ivGhost = a_data.ghostVect();

  Interval interv(0, m_nComp-1);

  for (int i = 0; i < m_qBar.size(); i++)
    {
      m_qBar[i] = new LevelData<EBCellFAB>(a_eblg.getDBL(), m_nComp, ivGhost, fact);
      m_qBar[i]->define(a_eblg.getDBL(), m_nComp, ivGhost, fact);
      a_data.copyTo(interv, *(m_qBar[i]), interv);

      if (m_doPIControl)
      {
        m_qdiffOld[i] = new LevelData<EBCellFAB>();
        m_qdiffNew[i] = new LevelData<EBCellFAB>();
        m_qdiffSum[i] = new LevelData<EBCellFAB>();

        m_qdiffOld[i]->define(a_eblg.getDBL(), m_nComp, ivGhost, fact);
        m_qdiffNew[i]->define(a_eblg.getDBL(), m_nComp, ivGhost, fact);
        m_qdiffSum[i]->define(a_eblg.getDBL(), m_nComp, ivGhost, fact);

        m_qBar_tOld[i] = new LevelData<EBCellFAB>();
        m_qBar_tOld[i]->define(a_eblg.getDBL(), m_nComp, ivGhost, fact);
      
        m_qdiffSum_tOld[i] = new LevelData<EBCellFAB>();
        m_qdiffSum_tOld[i]->define(a_eblg.getDBL(), m_nComp, ivGhost, fact);

        EBLevelDataOps::setToZero(*m_qdiffOld[i]);
        EBLevelDataOps::setToZero(*m_qdiffNew[i]);
        EBLevelDataOps::setToZero(*m_qdiffSum[i]);

        EBLevelDataOps::setToZero(*m_qBar_tOld[i]);
        EBLevelDataOps::setToZero(*m_qdiffSum_tOld[i]);
      }
    }
}
/*********/
void ChomboSFDInterface::
postInitialize(ChomboSFDInterface* a_finerLevel,
               EBCoarseAverage& a_ebCoarseAverage)
{
  doEBCoarseAverage(a_finerLevel, a_ebCoarseAverage);
}
/*********/
void ChomboSFDInterface::
doEBCoarseAverage(ChomboSFDInterface* a_finerLevel,
                  EBCoarseAverage& a_ebCoarseAverage)
{
  if (a_finerLevel != NULL)
  {
    Interval interv(0, m_nComp-1);

    for (int i = 0; i < m_qBar.size(); i++)
      {
        a_ebCoarseAverage.average(*(m_qBar[i]), *((a_finerLevel->m_qBar)[i]), interv);

        if (m_doPIControl)
        {
          a_ebCoarseAverage.average(*(m_qdiffOld[i]), *((a_finerLevel->m_qdiffOld)[i]), interv);
          a_ebCoarseAverage.average(*(m_qdiffNew[i]), *((a_finerLevel->m_qdiffNew)[i]), interv);
          a_ebCoarseAverage.average(*(m_qdiffSum[i]), *((a_finerLevel->m_qdiffSum)[i]), interv);

          a_ebCoarseAverage.average(*(m_qBar_tOld[i]), *((a_finerLevel->m_qBar_tOld)[i]), interv);
          a_ebCoarseAverage.average(*(m_qdiffSum_tOld[i]), *((a_finerLevel->m_qdiffSum_tOld)[i]), interv);          
        }
      }
  }
}
/*********/
void ChomboSFDInterface::
operator()(LevelData<EBCellFAB>& a_q, double a_dt, EBLevelGrid& a_eblg, double a_time, int a_curStep)
{
  CH_assert(m_isDefined);
  
  ChomboSFDOpFunctions functions;
  functions.setEBLG(&a_eblg);

  SFDOpFunctions<LevelData<EBCellFAB> >* cast_func = static_cast<SFDOpFunctions<LevelData<EBCellFAB> >* >(&functions);

  base_aXbY_ cast_aXbY = static_cast<base_aXbY_>(&ChomboSFDOpFunctions::aXbY);
  base_copier_ cast_copier = static_cast<base_copier_>(&ChomboSFDOpFunctions::copierFunc);
 
  m_SFDOp(a_q, m_qBar.stdVector(), a_dt, cast_aXbY, cast_copier, cast_func); 

  // for PI control
  // update qbarNew
  if (m_doPIControl)
  {
    Interval interv(0, m_nComp-1);
    for (int i = 0; i < m_qBar.size(); i++)
    {
      EBLevelDataOps::setToZero(*m_qdiffNew[i]);
      EBLevelDataOps::axby(*m_qdiffNew[i], a_q, *m_qBar[i], 1., -1.);

      // set source = -sum(integralCoef*integral(qdiff))
      //    integral(qdiff) from t1 to t2 =  a_dt*0.5*(qdiffNew + qdiffOld)
      EBLevelDataOps::incr(*m_qdiffSum[i], *m_qdiffNew[i], a_dt/2.);
      EBLevelDataOps::incr(*m_qdiffSum[i], *m_qdiffOld[i], a_dt/2.);
    }
  }
}
/*********/
void ChomboSFDInterface::
resetIntegrator(LevelData<EBCellFAB>& a_q, int a_iFilter, bool a_turnOff)
{
  if (m_doPIControl)
  {
//    m_integralCoef[a_iFilter] *= 1.1;
    EBLevelDataOps::setToZero(*m_qdiffOld[a_iFilter]);
    EBLevelDataOps::setToZero(*m_qdiffNew[a_iFilter]);
    EBLevelDataOps::setToZero(*m_qdiffSum[a_iFilter]);
    m_nResetSteps++; 
  }
}
/*********/
void ChomboSFDInterface::
resetIntegratorCoef(Real a_resetVal, int a_filterIndex)
{
  m_integralCoef[a_filterIndex] *= a_resetVal; 
}
/*********/
void ChomboSFDInterface::
resetState(const LevelData<EBCellFAB>& a_val, int a_iFilter)
{
  EBLevelDataOps::setToZero(*m_qBar[a_iFilter]);
  EBLevelDataOps::incr(*m_qBar[a_iFilter], a_val, 1.);

/*
  int nvar = a_q.nComp();
  const DisjointBoxLayout& grids = a_q.getBoxes();
  DataIterator dit = grids.dataIterator();
  // iterator over the grids on this processor
  for (dit.begin(); dit.ok(); ++dit)
  {
    EBCellFAB& q = a_q[dit()];
    EBCellFAB& qBar = (*m_qBar[a_iFilter])[dit()];
    EBCellFAB& qBar_tOld = (*m_qBar_tOld[a_iFilter])[dit()];
    EBCellFAB& qdiffSum = (*m_qdiffSum[a_iFilter])[dit()];
    EBCellFAB& qdiffSum_tOld = (*m_qdiffSum_tOld[a_iFilter])[dit()];

    const Box& box = grids.get(dit());
    const EBISBox& ebisBox = q.getEBISBox();
    IntVectSet ivs(box);
    for (VoFIterator vofit(ivs, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();

      Real tmp, ratio;
      for (int ivar = 0; ivar < nvar; ivar++)
      {  
        tmp = qBar(vof, ivar);
        ratio = (qBar(vof, ivar) - qBar_tOld(vof, ivar))/(qdiffSum(vof, ivar) - qdiffSum_tOld(vof, ivar));
        if (abs(ratio) > 1.) 
        {
          int sign = (ratio < 0) ? -1 : 1;
          ratio = 0.1*sign;
        }
        qBar(vof, ivar) -= 1.*qdiffSum(vof, ivar)*ratio;
//          qBar(vof,ivar) += 0.1;
        q(vof, ivar) = qBar(vof, ivar);
        qBar_tOld(vof, ivar) = tmp;
        qdiffSum_tOld(vof, ivar) = qdiffSum(vof, ivar);
      }
    }
  }
*/  
}
/*********/
void ChomboSFDInterface::
regrid(LevelData<EBCellFAB>& a_tempFAB,
       ChomboSFDInterface*   a_coarLevel,
       EBLevelGrid&          a_eblg,
       EBPWLFineInterp&      a_interp)
{
  CH_assert(m_isDefined);
 
  IntVect ivGhost = a_eblg.getGhost()*IntVect::Unit;
  EBCellFactory fact(a_eblg.getEBISL());

  for (int i = 0; i < m_qBar.size(); i++)
    {
      int nComp = m_qBar[i]->nComp();
      Interval interv(0, nComp-1);
      m_qBar[i]->copyTo(interv, a_tempFAB, interv);
      m_qBar[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);

      if (a_coarLevel != NULL)
        {
          a_interp.interpolate(*(m_qBar[i]), *((a_coarLevel->m_qBar)[i]), interv);
        }

      a_tempFAB.copyTo(interv, *(m_qBar[i]), interv);

      // regrid PI control stuff
      if (m_doPIControl)
      {
        m_qdiffOld[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffOld[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        if (a_coarLevel != NULL)
        {
          a_interp.interpolate(*(m_qdiffOld[i]), *((a_coarLevel->m_qdiffOld)[i]), interv);
        }
        a_tempFAB.copyTo(interv, *(m_qdiffOld[i]), interv);

        m_qdiffNew[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffNew[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        if (a_coarLevel != NULL)
        {
          a_interp.interpolate(*(m_qdiffNew[i]), *((a_coarLevel->m_qdiffNew)[i]), interv);
        }
        a_tempFAB.copyTo(interv, *(m_qdiffNew[i]), interv);

        m_qdiffSum[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffSum[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        if (a_coarLevel != NULL)
        {
          a_interp.interpolate(*(m_qdiffSum[i]), *((a_coarLevel->m_qdiffSum)[i]), interv);
        }
        a_tempFAB.copyTo(interv, *(m_qdiffSum[i]), interv);

        m_qBar_tOld[i]->copyTo(interv, a_tempFAB, interv);
        m_qBar_tOld[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        if (a_coarLevel != NULL)
        {
          a_interp.interpolate(*(m_qBar_tOld[i]), *((a_coarLevel->m_qBar_tOld)[i]), interv);
        }
        a_tempFAB.copyTo(interv, *(m_qBar_tOld[i]), interv);

        m_qdiffSum_tOld[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffSum_tOld[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);

        if (a_coarLevel != NULL)
        {
          a_interp.interpolate(*(m_qdiffSum_tOld[i]), *((a_coarLevel->m_qdiffSum_tOld)[i]), interv);
        }

        a_tempFAB.copyTo(interv, *(m_qdiffSum_tOld[i]), interv);
      }
    }
}
/*********/
void ChomboSFDInterface::
regridBaseData(LevelData<EBCellFAB>&                 a_tempFAB,
               EBLevelGrid&                          a_eblg)
{ 
  IntVect ivGhost = a_eblg.getGhost()*IntVect::Unit;
  EBCellFactory fact(a_eblg.getEBISL());

  for (int i = 0; i < m_qBar.size(); i++)
    {
      int nComp = m_qBar[i]->nComp();
      Interval interv(0, nComp-1);
      m_qBar[i]->copyTo(interv, a_tempFAB, interv);
      m_qBar[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);

      a_tempFAB.copyTo(interv, *(m_qBar[i]), interv);

      // regrid PI control stuff
      if (m_doPIControl)
      {
        m_qdiffOld[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffOld[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        a_tempFAB.copyTo(interv, *(m_qdiffOld[i]), interv);

        m_qdiffNew[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffNew[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        a_tempFAB.copyTo(interv, *(m_qdiffNew[i]), interv);

        m_qdiffSum[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffSum[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        a_tempFAB.copyTo(interv, *(m_qdiffSum[i]), interv);

        m_qBar_tOld[i]->copyTo(interv, a_tempFAB, interv);
        m_qBar_tOld[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        a_tempFAB.copyTo(interv, *(m_qBar_tOld[i]), interv);

        m_qdiffSum_tOld[i]->copyTo(interv, a_tempFAB, interv);
        m_qdiffSum_tOld[i]->define(a_eblg.getDBL(), nComp, ivGhost, fact);
        a_tempFAB.copyTo(interv, *(m_qdiffSum_tOld[i]), interv);
      }
    } 
}
/*********/
Vector<string> ChomboSFDInterface::
getDataNames() const
{
  Vector<string> retval;
  for (int i=0; i < m_qBar.size(); i++)
    {
      Vector<string> temp = getDataNames(i);
//      retval.insert(retval.end(), temp.begin(), temp.end());
      for (int comp = 0; comp < temp.size(); comp++)
        {
          retval.push_back(temp[comp]);
        }
    }

  return retval;
}
/*********/
Vector<string> ChomboSFDInterface::
getDataNames(int a_filterIndex) const
{
  Vector<string> retval;
  string index = SSTR(a_filterIndex);
  retval.push_back("qbar" + index + "0");
  retval.push_back("qbar" + index + "1");
  if (SpaceDim >= 2)
    {
      retval.push_back("qbar" + index + "2");
    }

  if (SpaceDim >= 3)
    {
      retval.push_back("qbar" + index + "3");
    }
  retval.push_back("qbar" + index + "4");

  retval.resize(m_nComp); // fix for INS

  return retval;
}
/*********/
int ChomboSFDInterface::
getnComp() const
{
  int retVal = m_nComp*m_qBar.size();
  return retVal;
}
/*********/
void ChomboSFDInterface::
appendDataToFAB(LevelData<FArrayBox>& a_result,
                int                   a_startIndex,
                const EBLevelGrid&          a_eblg) const
{
  CH_assert(m_isDefined);

  int nComp = 0;

  for (int i = 0; i < m_qBar.size(); i++ )
    {
      int startIndex = a_startIndex + nComp;
      appendDataToFAB(a_result, *(m_qBar[i]), startIndex, a_eblg);
      nComp += m_qBar[i]->nComp();
    }

}
/*********/
void ChomboSFDInterface::
appendDataToFAB(LevelData<FArrayBox>&       a_result,
                const LevelData<EBCellFAB>& a_data,
                int                         a_startIndex,
                const EBLevelGrid&                a_eblg) const
{
  CH_assert(m_isDefined)
 
  Real coveredVal = -1.2345678e-9;

  int nComp = a_data.nComp();

  for (DataIterator dit = a_eblg.getDBL().dataIterator(); dit.ok(); ++dit)
    {
      const EBISBox& ebisbox = a_eblg.getEBISL()[dit()];
      const Box& grid = a_eblg.getDBL().get(dit());

      const EBCellFAB& dataFab = a_data[dit()];

      FArrayBox& currentFab = a_result[dit()];

      // copy regular data
      currentFab.copy(dataFab.getSingleValuedFAB(), 0, a_startIndex, nComp);

      for (BoxIterator bit(grid); bit.ok(); ++bit)
        {
          const IntVect& iv = bit();
          // set special values for covered cells
          if (ebisbox.isCovered(iv))
            {
              for (int icomp = 0; icomp < nComp; icomp++)
                {
                  currentFab(iv, a_startIndex+icomp) = coveredVal;
                }
            }
        } 
    }
}
/*********/
const Vector<LevelData<EBCellFAB>* >* ChomboSFDInterface::
getData() const
{
  return &m_qBar;
}
/*********/
const LevelData<EBCellFAB>* ChomboSFDInterface::
getData(int a_filterIndex) const
{
  return m_qBar[a_filterIndex];
}
/*********/
const Vector<LevelData<EBCellFAB>* >* ChomboSFDInterface::
getIntegratedError() const
{
  if (m_doPIControl) 
  {
    return &m_qdiffSum;
  }
  else 
  {
    return NULL;
  }
}
/*********/
const Vector<LevelData<EBCellFAB>* >*  ChomboSFDInterface::
getQdiffNew() const
{
  if (m_doPIControl)
  {
    return &m_qdiffNew;
  }
  else
  {
    return NULL;
  }
}
/*********/
const Vector<LevelData<EBCellFAB>* >*  ChomboSFDInterface::
getQdiffOld() const
{
  if (m_doPIControl)
  {
    return &m_qdiffOld;
  }
  else
  {
    return NULL;
  }
}
/*********/
const LevelData<EBCellFAB>*  ChomboSFDInterface::
getQdiffNew(int a_filterIndex) const
{
  if (m_doPIControl)
  {
    return m_qdiffNew[a_filterIndex];
  }
  else
  {
    return NULL;
  }
}
/*********/
const LevelData<EBCellFAB>*  ChomboSFDInterface::
getQdiffOld(int a_filterIndex) const
{
  if (m_doPIControl)
  {
    return m_qdiffOld[a_filterIndex];
  }
  else
  {
    return NULL;
  }
}
/*********/
const LevelData<EBCellFAB>*  ChomboSFDInterface::
getIntegratedError(int a_filterIndex) const
{
  if (m_doPIControl)
  {
    return m_qdiffSum[a_filterIndex];
  }
  else
  {
    return NULL;
  }
}
/*********/
int ChomboSFDInterface::
getnFilters() const
{
  return m_qBar.size();
}
/*********/
bool ChomboSFDInterface::
convergedToSteadyState(const LevelData<EBCellFAB>& a_data, EBLevelGrid& a_eblg)
{
  EBCellFactory      cellFact(a_eblg.getEBISL());
  LevelData<EBCellFAB>  udiff(a_eblg.getDBL(), a_data.nComp(), 4*IntVect::Unit, cellFact);
  EBLevelDataOps::setToZero(udiff);
  EBLevelDataOps::incr(udiff, a_data, 1.0);

  // Sid: hardwiring it to 1st filter
  EBLevelDataOps::incr(udiff, *(m_qBar[0]),  -1.0);
  Real umax, umin, eps;
  Real dmax, dmin;
  int ivar;
  ParmParse pp;
  pp.get("convergence_metric", eps);
  pp.get("convergence_variable", ivar);
  EBLevelDataOps::getMaxMin(dmax, dmin, udiff, ivar);
//  EBLevelDataOps::getMaxMin(umax, umin, a_data, ivar);
  Real denom = 1;
/*
  if(Abs(umax - umin) > eps)
    {
      denom = Abs(umax-umin);
    }
*/
  Real maxdiff = Abs(dmax-dmin)/denom;
  pout() << "max difference in convergence variable = " << maxdiff << ", eps set to " << eps << endl;
  return (maxdiff < eps);
}
/*********/
void ChomboSFDInterface::
getNavierStokesSource(LevelData<EBCellFAB>& a_source, Real a_dt, bool a_updateOldToNew)
{
  EBLevelDataOps::setToZero(a_source);
  if (m_doPIControl)
  {
    CH_assert(a_source.nComp() == m_qdiffOld[0]->nComp());

    // set source = -sum(integralCoef*integral(qdiff))
    //    integral(qdiff) from t1 to t2 =  a_dt*0.5*(qdiffNew + qdiffOld)
    for (int i = 0; i < m_qBar.size(); i++)
    {
/*
      EBLevelDataOps::incr(a_source, *m_qdiffNew[i], -m_integralCoef[i]*a_dt/2.);
      EBLevelDataOps::incr(a_source, *m_qdiffOld[i], -m_integralCoef[i]*a_dt/2.);
*/
      EBLevelDataOps::incr(a_source, *m_qdiffSum[i], -m_integralCoef[i]);
    }

    if (a_updateOldToNew)
    {
      for (int i = 0; i < m_qBar.size(); i++)
      {
        EBLevelDataOps::setToZero(*m_qdiffOld[i]);
        EBLevelDataOps::incr(*m_qdiffOld[i], *m_qdiffNew[i], 1.);
      }
    }
  }
}
/*********/
void ChomboSFDInterface::
getNavierStokesSource(LevelData<EBCellFAB>& a_source, const int& a_startSrcComp, const int& a_startSFDComp, const int& a_nComp, Real a_dt, bool a_updateOldToNew)
{
  CH_assert(a_source.nComp() == 1);
  EBLevelDataOps::setToZero(a_source);
  if (m_doPIControl)
  {
    // set source = -sum(integralCoef*integral(qdiff))
    //    integral(qdiff) from t1 to t2 =  a_dt*0.5*(qdiffNew + qdiffOld)
    for (int i = 0; i < m_qBar.size(); i++)
    {
      for (DataIterator dit = a_source.dataIterator(); dit.ok(); ++dit)
      {
        DataIndex d = dit();
/*
        a_source[d].plus((*m_qdiffNew[i])[d], a_startSFDComp, a_startSrcComp, a_nComp);
        a_source[d].plus((*m_qdiffOld[i])[d], a_startSFDComp, a_startSrcComp, a_nComp);
        a_source[d].mult(-m_integralCoef[i]*a_dt/2.);
*/
        a_source[d].plus((*m_qdiffSum[i])[d], a_startSFDComp, a_startSrcComp, a_nComp);
        a_source[d].mult(-m_integralCoef[i]);
      }
    }

    if (a_updateOldToNew)
    {
      for (int i = 0; i < m_qBar.size(); i++)
      {
        EBLevelDataOps::setToZero(*m_qdiffOld[i]);
        EBLevelDataOps::incr(*m_qdiffOld[i], *m_qdiffNew[i], 1.);
      }
    }
  }
}
/*********/
#ifdef CH_USE_HDF5

/*********/
void ChomboSFDInterface::
writeCheckpointLevel(HDF5Handle& a_handle) const
{
  for (int i = 0; i < m_qBar.size(); i++)
    {
      writeCheckpointLevel(a_handle, i);
    }
}
/*********/
void ChomboSFDInterface::
writeCheckpointLevel(HDF5Handle& a_handle, int a_filterIndex) const
{
  string str1 = "qBar"+SSTR(a_filterIndex);
  write(a_handle, *(m_qBar[a_filterIndex]), str1);

  if (m_doPIControl)
  {
    string str2 = "qdiffOld"+SSTR(a_filterIndex);
    write(a_handle, *(m_qdiffOld[a_filterIndex]), str2);

    string str3 = "qdiffNew"+SSTR(a_filterIndex);
    write(a_handle, *(m_qdiffNew[a_filterIndex]), str3);

    string str4 = "qdiffSum"+SSTR(a_filterIndex);
    write(a_handle, *(m_qdiffSum[a_filterIndex]), str4);
  }
}
/*********/
void ChomboSFDInterface::
readCheckpointLevel(HDF5Handle& a_handle, EBLevelGrid& a_eblg)
{

  for (int i = 0; i < m_qBar.size(); i++)
    {
      readCheckpointLevel(a_handle, i, a_eblg);
    }
}
/*********/
void ChomboSFDInterface::
readCheckpointLevel(HDF5Handle& a_handle, int a_filterIndex, EBLevelGrid& a_eblg) 
{
  CH_assert(m_isDefined);

  EBCellFactory factoryNew(a_eblg.getEBISL());
  const IntVect ivGhost = a_eblg.getGhost()*IntVect::Unit;
 
  m_qBar[a_filterIndex] = new LevelData<EBCellFAB>(a_eblg.getDBL(), m_nComp, ivGhost, factoryNew);

  m_qBar[a_filterIndex]->define(a_eblg.getDBL(), m_nComp, ivGhost, factoryNew);

  if (m_doPIControl)
  {
    m_qdiffOld[a_filterIndex] = new LevelData<EBCellFAB>();
    m_qdiffNew[a_filterIndex] = new LevelData<EBCellFAB>();
    m_qdiffSum[a_filterIndex] = new LevelData<EBCellFAB>();

    m_qdiffOld[a_filterIndex]->define(a_eblg.getDBL(), m_nComp, ivGhost, factoryNew);
    m_qdiffNew[a_filterIndex]->define(a_eblg.getDBL(), m_nComp, ivGhost, factoryNew);
    m_qdiffSum[a_filterIndex]->define(a_eblg.getDBL(), m_nComp, ivGhost, factoryNew);
  }

  string str1 = "qBar"+SSTR(a_filterIndex);

  // the false says to not redefine data
  int dataStatus = read<EBCellFAB>(a_handle, *(m_qBar[a_filterIndex]), str1, a_eblg.getDBL(), Interval(), false);

  if ((dataStatus != 0))
    {
      MayDay::Error("file does not contain qBar data");
    }

  if (m_doPIControl)
  {
    string str2 = "qdiffOld"+SSTR(a_filterIndex);
    string str3 = "qdiffNew"+SSTR(a_filterIndex);
    string str4 = "qdiffSum"+SSTR(a_filterIndex);

    dataStatus = read<EBCellFAB>(a_handle, *(m_qdiffOld[a_filterIndex]), str2, a_eblg.getDBL(), Interval(), false);
    if ((dataStatus != 0))
    {
      MayDay::Error("file does not contain qdiffOld data");
    }

    dataStatus = read<EBCellFAB>(a_handle, *(m_qdiffNew[a_filterIndex]), str3, a_eblg.getDBL(), Interval(), false);
    if ((dataStatus != 0))
    {
      MayDay::Error("file does not contain qdiffNew data");
    }

    dataStatus = read<EBCellFAB>(a_handle, *(m_qdiffSum[a_filterIndex]), str4, a_eblg.getDBL(), Interval(), false);
    if ((dataStatus != 0))
    {
      MayDay::Error("file does not contain qdiffNew data");
    }
  }
}
/*********/
void ChomboSFDInterface::
writePlotHeader(HDF5HeaderData& a_header, int a_startIndex) const
{
  Vector<string> names = getDataNames();
  a_header.m_int["num_components"] = a_startIndex + names.size();

  char compStr[30];
  for (int comp = a_startIndex; comp < (a_startIndex+names.size()); ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      a_header.m_string[compStr] = names[comp-a_startIndex];
    }
}
/*********/
void ChomboSFDInterface::
writePlotHeader(HDF5HeaderData& a_header, int a_comp, int a_startIndex) const
{
  Vector<string> names = getDataNames(a_comp);
  a_header.m_int["num_components"] = a_startIndex + names.size();

  char compStr[30];
  for (int comp = a_startIndex; comp < (a_startIndex+names.size()); ++comp)
    {
      sprintf(compStr,"component_%d",comp);
      a_header.m_string[compStr] = names[comp-a_startIndex];
    }
}
/*********/
#endif
