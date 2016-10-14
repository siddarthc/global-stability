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
}
/*********/
ChomboSFDInterface::~ChomboSFDInterface()
{

  for (int i = 0; i< m_qBar.size(); i++)
    {
      if (m_qBar[i] != NULL) delete m_qBar[i];
    }

}
/*********/
void ChomboSFDInterface::
define(int    a_nFilters,
       double a_smallestFilter,
       double a_largestFilter,
       double a_controlCoef,
       int    a_nComp)
{
  m_qBar.resize(a_nFilters);
  m_SFDOp.define(a_nFilters, a_smallestFilter, a_largestFilter, a_controlCoef);
  m_nComp = a_nComp;
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
      }
  }
}
/*********/
void ChomboSFDInterface::
operator()(LevelData<EBCellFAB>& a_q, double a_dt, EBLevelGrid& a_eblg)
{
  CH_assert(m_isDefined);
  
  ChomboSFDOpFunctions functions;
  functions.setEBLG(&a_eblg);

  SFDOpFunctions<LevelData<EBCellFAB> >* cast_func = static_cast<SFDOpFunctions<LevelData<EBCellFAB> >* >(&functions);

  base_aXbY_ cast_aXbY = static_cast<base_aXbY_>(&ChomboSFDOpFunctions::aXbY);
  base_copier_ cast_copier = static_cast<base_copier_>(&ChomboSFDOpFunctions::copierFunc);
 
  m_SFDOp(a_q, m_qBar.stdVector(), a_dt, cast_aXbY, cast_copier, cast_func); 
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
  return m_nComp*m_qBar.size();
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
Vector<LevelData<EBCellFAB>* >* ChomboSFDInterface::
getData()
{
  return &m_qBar;
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
  string str = "qBar"+SSTR(a_filterIndex);
  write(a_handle, *(m_qBar[a_filterIndex]), str);
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

  string str = "qBar"+SSTR(a_filterIndex);

  // the false says to not redefine data
  int dataStatus = read<EBCellFAB>(a_handle, *(m_qBar[a_filterIndex]), str, a_eblg.getDBL(), Interval(), false);

  if ((dataStatus != 0))
    {
      MayDay::Error("file does not contain qBar data");
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
