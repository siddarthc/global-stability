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
  //
  static void rollChomboDataToEpetraVec(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                        Epetra_Vector*                       a_EpetraVec,
                                        bool                                 a_incOverlapData = false);


  // function to put Epetra_Vector's data to Chombo's data Level by Level
  //   Assumes that the data in Epetra_Vector of the structure: 
  //     finestLevelData followed by its coarser level so on to coarsest level
  //
  // if a_isOverlapDataInc is a flag to indicate if Epetra_Vector has the data at the Coarse-Fine overlap regions incluced in it
  //    if true, Epetra_Vector has the overlap data included
  //
  static void unrollEpetraVecToChomboData(Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                          const Epetra_Vector*                 a_EpetraVec,
                                          bool                           a_isOverlapDataInc = false);


  // get the number of Elements of a_ChomboData on the current proc
  // if a_incOverlapData is true, data at the Coarse-Fine overlap is also counted
  //
  static int getnElementsOnThisProc(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                    bool a_incOverlapData = false);


  // get the Epetra_Map corresponding to the given a_ChomboData and a_commPtr
  // if a_incOverlapData is true, data at the Coarse-Fine overlap is also counted
  static Epetra_Map getEpetraMap(const Vector<LevelData<EBCellFAB>* >& a_ChomboData,
                                 const Epetra_Comm*                    a_commPtr,
                                 bool a_incOverlapData = false);


};
#endif
