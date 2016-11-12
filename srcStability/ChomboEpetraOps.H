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

///
/**
  class to manage CHOMBO Epetra Operations
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
                                        const Vector<Real>&                  a_refRatio = Vector<Real>());


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
                                          const Vector<Real>&            a_refRatio = Vector<Real>(),
                                          double                         a_overlapVal = 0.);


  // get the number of Elements of a_ChomboData on the current proc
  // if a_incOverlapData is true, data at the Coarse-Fine overlap is also counted
  // a_refRatio is the refinement ratio between the levels. It is needed only if a_incOverlapData is false
  //
  static int getnElementsOnThisProc(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                    bool                                  a_incOverlapData = true,
                                    const Vector<Real>&                   a_refRatio = Vector<Real>());


  // get the Epetra_Map corresponding to the given a_ChomboData and a_commPtr
  // if a_incOverlapData is true, data at the Coarse-Fine overlap is also counted
  // a_refRatio is the refinement ratio between the levels. It is needed only if a_incOverlapData is false
  //
  static Epetra_Map getEpetraMap(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                 const Epetra_Comm*                    a_commPtr,
                                 bool                                  a_incOverlapData = true,
                                 const Vector<Real>&                   a_refRatio = Vector<Real>());

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

  static void getVolWeights(Epetra_Vector*                       a_weights,
                            double&                              a_domainVolume,
                            const Vector<LevelData<EBCellFAB>* > a_ChomboData, 
                            const double&                        a_coarsestDx,
                            const Vector<Real>&                  a_refRatio,
                            bool                                 a_incOverlapData = true);

  // compute L2 norm of the Epetra_Vector
  //  just computes the regular unweighted norm
  static int computeL2Norm(double&              a_norm,
                           const Epetra_Vector& a_mv);
/*
  // compute volume weighted L2 norm of the Epetra_Vector
  // a_n
  static int computeVolWeightedL2Norm(double&              a_norm,
                                      const Epetra_Vector& a_mv,
                                      double               a_coarsestDx,
                                      const Vector<Real>&  a_refRatio,
                                      const Vector<LevelData<EBCellFAB> >);
*/
};
#endif
