/*
 *
 *
 *
 *
 */

#include "AnasaziEpetraSolverAdapter.H"

namespace Anasazi {

		/* --- Construction/Destruction ---*/

  EpetraMultiVecSolverExt::EpetraMultiVecSolverExt(const Epetra_BlockMap& Map_in, double * array, const int numvecs, const RCP<TrilinosSolverInterface>& a_solverInterface, const int stride)
    : Epetra_MultiVector(Epetra_DataAccess::Copy, Map_in, array, stride, numvecs) 
  {
    m_solverInterface = a_solverInterface;
  }

  EpetraMultiVecSolverExt::EpetraMultiVecSolverExt(const Epetra_BlockMap& Map_in, const int numvecs, const RCP<TrilinosSolverInterface>& a_solverInterface)
    : Epetra_MultiVector(Map_in, numvecs) 
  {
    m_solverInterface = a_solverInterface;
  }
   
  EpetraMultiVecSolverExt::EpetraMultiVecSolverExt(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index, const RCP<TrilinosSolverInterface>& a_solverInterface)
    : Epetra_MultiVector(CV, P_vec, &(const_cast<std::vector<int> &>(index))[0], index.size())
  {
    m_solverInterface = a_solverInterface;
  }
  
  
  EpetraMultiVecSolverExt::EpetraMultiVecSolverExt(const Epetra_MultiVector& P_vec, const RCP<TrilinosSolverInterface>& a_solverInterface)
    : Epetra_MultiVector(P_vec) 
  {
    m_solverInterface = a_solverInterface;
  }

		/* --- Member functions inherited from Anasazi::MultiVec --- */

  
  MultiVec<double>* EpetraMultiVecSolverExt::Clone ( const int numvecs ) const
  {
    EpetraMultiVecSolverExt * ptr_apv = new EpetraMultiVecSolverExt(Map(), numvecs, m_solverInterface);
    return ptr_apv; // safe upcast.
  }

  MultiVec<double>* EpetraMultiVecSolverExt::CloneCopy() const
  {
    EpetraMultiVecSolverExt *ptr_apv = new EpetraMultiVecSolverExt(*this, m_solverInterface);
    return ptr_apv; // safe upcast
  }
  
  
  MultiVec<double>* EpetraMultiVecSolverExt::CloneCopy ( const std::vector<int>& index ) const
  {
    EpetraMultiVecSolverExt * ptr_apv = new EpetraMultiVecSolverExt(Epetra_DataAccess::Copy, *this, index, m_solverInterface);
    return ptr_apv; // safe upcast.
  }
  
  
  MultiVec<double>* EpetraMultiVecSolverExt::CloneViewNonConst ( const std::vector<int>& index ) 
  {
    EpetraMultiVecSolverExt * ptr_apv = new EpetraMultiVecSolverExt(Epetra_DataAccess::View, *this, index, m_solverInterface);
    return ptr_apv; // safe upcast.
  }
  
  const MultiVec<double>* EpetraMultiVecSolverExt::CloneView ( const std::vector<int>& index ) const
  {
    EpetraMultiVecSolverExt * ptr_apv = new EpetraMultiVecSolverExt(Epetra_DataAccess::View, *this, index, m_solverInterface);
    return ptr_apv; // safe upcast.
  }


    void EpetraMultiVecSolverExt::SetBlock( const MultiVec<double>& A, const std::vector<int>& index ) 
  {
    // this should be revisited to e
    EpetraMultiVecSolverExt temp_vec(Epetra_DataAccess::View, *this, index, m_solverInterface);

    int numvecs = index.size();
    if ( A.GetNumberVecs() != numvecs ) {
      std::vector<int> index2( numvecs );
      for(int i=0; i<numvecs; i++)
        index2[i] = i;
      EpetraMultiVecSolverExt *tmp_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(A)); 
      TEUCHOS_TEST_FOR_EXCEPTION( tmp_vec==NULL, std::invalid_argument, "Anasazi::EpetraMultiVecSolverExt::SetBlocks() cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");
      EpetraMultiVecSolverExt A_vec(Epetra_DataAccess::View, *tmp_vec, index2, m_solverInterface);
      temp_vec.MvAddMv( 1.0, A_vec, 0.0, A_vec );
    }
    else {
      temp_vec.MvAddMv( 1.0, A, 0.0, A );
    }
  }

  
  void EpetraMultiVecSolverExt::MvTimesMatAddMv ( double alpha, const MultiVec<double>& A, 
      const Teuchos::SerialDenseMatrix<int,double>& B, double beta ) 
  {
    Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
    Epetra_MultiVector B_Pvec(Epetra_DataAccess::View, LocalMap, B.values(), B.stride(), B.numCols());
    
    EpetraMultiVecSolverExt *A_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVecSolverExt::SetBlocks() cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");
    
    TEUCHOS_TEST_FOR_EXCEPTION( 
        Multiply( 'N', 'N', alpha, *A_vec, B_Pvec, beta ) != 0,
        EpetraMultiVecFailure, "Anasazi::EpetraMultiVecSolverExt::MvTimesMatAddMv() call to Epetra_MultiVec::Multiply() returned a nonzero value.");
  }


  void EpetraMultiVecSolverExt::MvAddMv ( double alpha , const MultiVec<double>& A, 
                                 double beta, const MultiVec<double>& B) 
  {
    EpetraMultiVecSolverExt *A_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVecSolverExt::MvAddMv() cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");
    EpetraMultiVecSolverExt *B_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(B)); 
    TEUCHOS_TEST_FOR_EXCEPTION( B_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVecSolverExt::MvAddMv() cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");
    
    TEUCHOS_TEST_FOR_EXCEPTION( 
        Update( alpha, *A_vec, beta, *B_vec, 0.0 ) != 0,
        EpetraMultiVecFailure, "Anasazi::EpetraMultiVecSolverExt::MvAddMv() call to Epetra_MultiVec::Update() returned a nonzero value.");
  }

  
  void EpetraMultiVecSolverExt::MvTransMv ( double alpha, const MultiVec<double>& A,
                                   Teuchos::SerialDenseMatrix<int,double>& B
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                                   , ConjType conj
#endif
                                  ) const
  {    
    EpetraMultiVecSolverExt *A_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(A));
    
    if (A_vec) {
      Epetra_LocalMap LocalMap(B.numRows(), 0, Map().Comm());
      Epetra_MultiVector B_Pvec(Epetra_DataAccess::View, LocalMap, B.values(), B.stride(), B.numCols());
      
    TEUCHOS_TEST_FOR_EXCEPTION( 
        B_Pvec.Multiply( 'T', 'N', alpha, *A_vec, *this, 0.0 ) != 0,
        EpetraMultiVecFailure, "Anasazi::EpetraMultiVecSolverExt::MvTransMv() call to Epetra_MultiVec::Multiply() returned a nonzero value.");
    }
  }


  void EpetraMultiVecSolverExt::MvDot ( const MultiVec<double>& A, std::vector<double> & b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
                               , ConjType conj
#endif
                             ) const
  {
    EpetraMultiVecSolverExt *A_vec = dynamic_cast<EpetraMultiVecSolverExt *>(&const_cast<MultiVec<double> &>(A)); 
    TEUCHOS_TEST_FOR_EXCEPTION( A_vec==NULL,  std::invalid_argument, "Anasazi::EpetraMultiVecSolverExt::MvDot() cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");

    if (( (int)b.size() >= A_vec->NumVectors() ) ) {
      TEUCHOS_TEST_FOR_EXCEPTION( 
          m_solverInterface->computeDotProd( *this, *A_vec, &b[0] ) != 0,
          EpetraMultiVecFailure, "Anasazi::EpetraMultiVecSolverExt::MvDot() call to Epetra_MultiVec::Dot() returned a nonzero value.");
    }
  }


  void EpetraMultiVecSolverExt::MvNorm ( std::vector<double> & normvec) const
  {
    if (((int)normvec.size() >= GetNumberVecs()) ) {
        TEUCHOS_TEST_FOR_EXCEPTION( m_solverInterface->computeL2Norm(*this, &normvec[0])!=0, EpetraMultiVecFailure, "Anasazi::EpetraMultiVecSolverExt::MvNorm call to Epetra_MultiVector::Norm2() returned a nonzero value.");
      }
  }

  void EpetraMultiVecSolverExt::MvScale ( const std::vector<double>& alpha )
  {
    // Check to make sure the vector is as long as the multivector has columns.
    int numvecs = this->NumVectors();
    TEUCHOS_TEST_FOR_EXCEPTION( (int)alpha.size() != numvecs, std::invalid_argument, 
        "Anasazi::EpetraMultiVecSolverExt::MvScale() alpha argument size was inconsistent with number of vectors in mv.");
    
    std::vector<int> tmp_index( 1, 0 );
    for (int i=0; i<numvecs; i++) {
      Epetra_MultiVector temp_vec(Epetra_DataAccess::View, *this, &tmp_index[0], 1);
      TEUCHOS_TEST_FOR_EXCEPTION( 
          temp_vec.Scale( alpha[i] ) != 0,
          EpetraMultiVecFailure, "Anasazi::EpetraMultiVecSolverExt::MvScale() call to Epetra_MultiVec::Scale() returned a nonzero value.");
      tmp_index[0]++;
    }
  }

} // end namespace Anasazi
