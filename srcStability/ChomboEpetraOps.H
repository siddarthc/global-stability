/*
 *
 *
 *
 *
 */

#ifndef _Chombo_Epetra_Ops_H_
#define _Chombo_Epetra_Ops_H_

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "LevelData.H"
#include "Vector.H"
#include "EBCellFAB.H"
#include "EBLevelGrid.H"

///
/**
  class to manage CHOMBO Epetra Operations
  CAUTION: These methods work only if Epetra uses same Map as CHOMBO
 */
class ChomboEpetraOps
{
public:
  ChomboEpetraOps()
  {
  }

  ~ChomboEpetraOps()
  {
  }

  // function to put Chombo's data on multiple levels into Epetra_Vector
  //   if a_ChomboData has 3 Comp (a,b,c), the structure of the EpetraVec is 
  //        : [(a_1, b_1 , c_1 , a_2, b_2, c_2, ... , a_n, b_n, c_n) 
  //           where 1,2,3,...,n are data indices
  //
  // On Multiple levels, the data is written level by level
  // Finest level data is written first, followed by its coarser level, so on to the coarsest level
  //
  // if a_incOverlapData is true, the data at the Coarse-Fine overlap regions is also copied
  // a_refRatio is the refinement ratio between the levels. It is needed only if a_incOverlapData is false
  //
  static void rollChomboDataToEpetraVec(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                        Epetra_Vector*                       a_EpetraVec,
                                        bool                                 a_incOverlapData = true,
                                        const Vector<int>&                   a_refRatio = Vector<int>());


  // function to put Epetra_Vector's data to Chombo's data Level by Level
  //   Assumes that the data in Epetra_Vector of the structure: 
  //     finestLevelData followed by its coarser level so on to coarsest level
  //
  // if a_isOverlapDataInc is a flag to indicate if Epetra_Vector has the data at the Coarse-Fine overlap regions incluced in it
  //    if true, Epetra_Vector has the overlap data included
  //  a_refRatio is the refinement ratio between the levels. It is needed only if a_overlapDataInc is false
  //  a_overlapVal is the value that is set at the region where coarse level overlaps with finer level. 
  //      This is required only if a_isOverlapDataInc = false
  //
  static void unrollEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                          const Epetra_Vector*           a_EpetraVec,
                                          bool                           a_isOverlapDataInc = true,
                                          const Vector<int>&             a_refRatio = Vector<int>(),
                                          double                         a_overlapVal = 0.);


  // makes a_v = scaleV*a_v + scaleCD*a_ChomboData
  // a_startCompV is the starting component of a_v
  // a_startCompCD is the starting comp of a_ChomboData
  // a_nComp is the number of components to be added
  // a_totComp is the total no. of components in the a_v (not the length.)
  // 	For example, if a_v contains velx,vely,velz,pres then totComp = 4
  // 		and if only the value corresponding to pressure needs to be added, nComp = 1 
  // 			and startCompV is 3
  static void addChomboDataToEpetraVec(Epetra_Vector*                        a_v,
                                       const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                       const double&                         a_scaleV,
                                       const double&                         a_scaleCD,
                                       const int&                            a_startCompV,
                                       const int&                            a_startCompCD,
                                       const int&                            a_nComp,
                                       const int&                            a_totComp,
                                       bool                                  a_incOverlapData = true,
                                       const Vector<int>&                    a_refRatio = Vector<int>());

  // get the number of Elements of a_ChomboData on the current proc
  // if a_incOverlapData is true, data at the Coarse-Fine overlap is also counted
  // a_refRatio is the refinement ratio between the levels. It is needed only if a_incOverlapData is false
  //
  static int getnElementsOnThisProc(const Vector<DisjointBoxLayout>& a_ChomboDBL,
                                    const Vector<EBLevelGrid>&       a_ChomboEBLG,
                                    const int&                       a_nComp,
                                    bool                             a_incOverlapData = true,
                                    const Vector<int>&               a_refRatio = Vector<int>());


  // get the Epetra_Map corresponding to the given a_ChomboData and a_commPtr
  // if a_incOverlapData is true, data at the Coarse-Fine overlap is also counted
  // a_refRatio is the refinement ratio between the levels. It is needed only if a_incOverlapData is false
  //
  static Epetra_Map getEpetraMap(const Vector<DisjointBoxLayout>& a_ChomboDBL,
                                 const Vector<EBLevelGrid>&       a_ChomboEBLG,
                                 const int&                       a_nComp,
                                 const Epetra_Comm*               a_commPtr,
                                 bool                             a_incOverlapData = true,
                                 const Vector<int>&               a_refRatio = Vector<int>());

  // get the weights for each component of the Epetra_Vector to compute weighted L2 norm and weighted inner product.
  // 
  // check the definition of weighted norm in Epetra_MultiVector doxygen
  // weight = sqrt(domain volume/(volfrac*cellvol*Epetra_Vector->globallength()) )
  //
  // a_volWeights is the vector of volfrac*cellArea of all the cells 
  // a_domainVolume is the total volume of the domain
  // a_chomboData is the data which defines the Epetra_Vector's structure
  // a_coarsestDx is the coarsest level's dx
  // a_refRatio is the vector of refinement ratios between levels
  // a_incOverlapData determined whether the coarse-fine overlap regions are included

  static void getVolWeights(Epetra_Vector*                   a_weights,
                            double&                          a_domainVolume,
                            const Vector<DisjointBoxLayout>& a_ChomboDBL,
                            const Vector<EBLevelGrid>&       a_ChomboEBLG,
                            const int&                       a_nComp, 
                            const double&                    a_coarsestDx,
                            const Vector<int>&               a_refRatio,
                            bool                             a_incOverlapData = true);

  // compute L2 norm of the Epetra_Vector
  //  just computes the regular unweighted norm
  static int computeL2Norm(double&              a_norm,
                           const Epetra_Vector& a_v);

  // compute volume weighted L2 norm of the Epetra_Vector
  // a_norm is the computed norm
  // a_v is the vector whose norm is computed
  // a_weights is the vector of weights
  static int computeWeightedL2Norm(double&              a_norm,
                                   const Epetra_Vector& a_v,
                                   const Epetra_Vector& a_weights);


  // compute the dot product between 2 vectors
  // a_v1 and a_v2 are the two vectors whose dot product is computed
  static int computeDotProduct(double&              a_result,
                               const Epetra_Vector& a_v1,
                               const Epetra_Vector& a_v2);

  // compute weights dot product of a_v1 and a_v2: 
  // 		a_result = <a_v1|a_weights|a_v2>
  // if a_unscale weights is true, it assumes weights are:
  // 		(V/(n*dv))^0.5, where V is domain volume, n is Global vector length and dv is cell volume
  // 		weights are unscaled and set to the cell volume dv
  // if a_unscale is false, it assumes incoming weights are correct
  static int computeWeightedDotProduct(double& a_result,
                                       const Epetra_Vector& a_v1,
                                       const Epetra_Vector& a_v2,
                                       const Epetra_Vector& a_weights,
                                       bool a_unscaleWeights = true,
                                       double a_domainVolume = 0.);
};
#endif
