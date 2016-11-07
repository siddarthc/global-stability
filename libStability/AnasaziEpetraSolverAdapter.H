/*
 *
 *
 *
 *
 */

#ifndef _ANASAZI_EPETRA_SOLVER_ADAPTER_H_
#define _ANASAZI_EPETRA_SOLVER_ADAPTER_H_

#include "AnasaziEpetraAdapter.hpp"
#include "TrilinosSolverInterface.H"
#include <Teuchos_RCPDecl.hpp>

using Teuchos::RCP;

namespace Anasazi {

  // Extension to the Anasazi::MultiVec that uses Epetra_MultiVec
  // Just changing the norm and inner product definitions
  // Everything is just a copy of the Anasazi::EpetraMultiVec

  class EpetraMultiVecSolverExt : public MultiVec<double>, public Epetra_MultiVector, public EpetraMultiVecAccessor 
  {
  public:
    		/* --- Constructors/Destructors --- */

    // Basic constructor
    EpetraMultiVecSolverExt(const Epetra_BlockMap& Map_in, const int numvecs, const RCP<TrilinosSolverInterface>& a_solverInterface);

    // Copy Constructor
    EpetraMultiVecSolverExt(const Epetra_MultiVector& P_vec, const RCP<TrilinosSolverInterface>& a_solverInterface);

    EpetraMultiVecSolverExt(const Epetra_BlockMap& Map_in, double* array, const int numvecs, const RCP<TrilinosSolverInterface>& a_solverInterface, const int strude = 0);

    EpetraMultiVecSolverExt(Epetra_DataAccess CV, const Epetra_MultiVector& P_vec, const std::vector<int>& index, const RCP<TrilinosSolverInterface>& a_solverInterface);

    virtual ~EpetraMultiVecSolverExt() {};

    // Creation methods
      
    MultiVec<double> * Clone ( const int numvecs ) const;

    // create new and copy contents of this into new (deep copy)
    // return pointer to the EpetraMultiVecSolverExt

    MultiVec<double> * CloneCopy () const;

    MultiVec<double> * CloneCopy ( const std::vector<int> & index) const;

    MultiVec<double> * CloneViewNonConst ( const std::vector<int>& index );

    const MultiVec<double> * CloneView ( const std::vector<int>& index ) const;

    ptrdiff_t GetGlobalLength () const 
    {
       if ( Map().GlobalIndicesLongLong() )
          return static_cast<ptrdiff_t>( GlobalLength64() );
       else 
          return static_cast<ptrdiff_t>( GlobalLength() );
    }

    int GetNumberVecs () const { return NumVectors(); }
 
    		/* --- update Methods: --- */
    //
    void MvTimesMatAddMv ( double alpha, const MultiVec<double>& A, 
                           const Teuchos::SerialDenseMatrix<int,double>& B, 
                           double beta );

    void MvAddMv ( double alpha, const MultiVec<double>& A, 
                   double beta, const MultiVec<double>& B);

    void MvTransMv ( double alpha, const MultiVec<double>& A, Teuchos::SerialDenseMatrix<int,double>& B 
#ifdef HAVE_ANASAZI_EXPERIMENTAL
        , ConjType conj = Anasazi::CONJ
#endif
        ) const;

    // Depends on the user defined inner-product definition
    void MvDot ( const MultiVec<double>& A, std::vector<double> &b
#ifdef HAVE_ANASAZI_EXPERIMENTAL
        , ConjType conj = Anasazi::CONJ
#endif
        ) const;

    void MvScale ( double alpha ) { 
      TEUCHOS_TEST_FOR_EXCEPTION( this->Scale( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVecSolverExt::MvScale call to Epetra_MultiVector::Scale() returned a nonzero value.");
    }

    void MvScale ( const std::vector<double>& alpha );

    // Depends on the user defined norm definition
    void MvNorm ( std::vector<double> & normvec ) const;

    		/* --- Initialization methods --- */

    void SetBlock ( const MultiVec<double>& A, const std::vector<int>& index );

    void MvRandom() { 
      TEUCHOS_TEST_FOR_EXCEPTION( this->Random()!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVecSolverExt::MvRandom call to Epetra_MultiVector::Random() returned a nonzero value.");
    };

    void MvInit ( double alpha ) { 
      TEUCHOS_TEST_FOR_EXCEPTION( this->PutScalar( alpha )!=0, EpetraMultiVecFailure,
          "Anasazi::EpetraMultiVecSolverExt::MvInit call to Epetra_MultiVector::PutScalar() returned a nonzero value.");
    };

    	/* --- Accessor Methods (inherited from EpetraMultiVecAccessor) --- */

    Epetra_MultiVector* GetEpetraMultiVec() { return this; };

    const Epetra_MultiVector* GetEpetraMultiVec() const { return this; };

    		/* --- Print Method --- */

    void MvPrint( std::ostream& os ) const { os << *this << std::endl; };

  protected:

    RCP<TrilinosSolverInterface> m_solverInterface;

  private:

  }; // end EpetraMultiVecSolverExt

} // end namespace Anasazi

#endif
