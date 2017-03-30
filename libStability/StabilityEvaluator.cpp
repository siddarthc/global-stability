/*
 *
 *
 *
 *
 */

#include "StabilityEvaluator.H"
#include "AnasaziEpetraSolverAdapter.H"

#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

/*********/
StabilityEvaluator::
StabilityEvaluator(double                              a_eps, 
                   double                              a_integrationTime,
                   string                              a_baseflowFile,
                   RCP<TrilinosSolverInterfaceFactory> a_solverInterfaceFact,
                   MPI_Comm*                           a_commPtr)
{
  define(a_eps, a_integrationTime, a_baseflowFile, a_solverInterfaceFact, a_commPtr);
}

/*********/
void StabilityEvaluator:: 
define(double                              a_eps, 
       double                              a_integrationTime,
       string                              a_baseflowFile,
       RCP<TrilinosSolverInterfaceFactory> a_solverInterfaceFact,
       MPI_Comm*                           a_commPtr)
{
  m_eps             = a_eps;
  m_integrationTime = a_integrationTime;
  m_baseflowFile    = a_baseflowFile;
  m_solverInterface = a_solverInterfaceFact->create();

  if (a_commPtr == NULL)
  {
    Epetra_SerialComm* commPtr = new Epetra_SerialComm ();
    m_commPtr = (Epetra_Comm*)commPtr;        
    m_isParallelRun = false;
  }
  else
  {
    Epetra_MpiComm* commPtr = new Epetra_MpiComm (*a_commPtr);
    m_commPtr = (Epetra_Comm*)commPtr;
    m_isParallelRun = true;
  }

  m_solverInterface->setEps(m_eps);
  m_solverInterface->setIntegrationTime(m_integrationTime);
  m_solverInterface->setBaseflow(m_baseflowFile, m_commPtr);

  m_isDefined       = true;
}

/*********/
StabilityEvaluator::
~StabilityEvaluator()
{
  delete m_commPtr;
}
/*********/
int StabilityEvaluator::
computeDominantModes(double      a_tol,
                     int         a_nev,
                     int         a_numBlocks,
                     int         a_blockSize,
                     int         a_maxRestarts,
                     string      a_which,
                     bool        a_isOpSymmetric,
                     bool        a_verbose,
                     vector<int> a_plotEVComps)
{
  int verbosity = Anasazi::Errors + Anasazi::Warnings + Anasazi::FinalSummary;

  Teuchos::ParameterList EigenPL;
  EigenPL.set( "Verbosity", verbosity);
  EigenPL.set( "Which", a_which);
  EigenPL.set( "Block Size", a_blockSize );
  EigenPL.set( "Num Blocks", a_numBlocks );
  EigenPL.set( "Maximum Restarts", a_maxRestarts );
  EigenPL.set( "Convergence Tolerance", a_tol );

  // just use the same distribution of elements as the solver to avoid comm overhead. Not sure if it's the best way
//  int locElements = m_solverInterface->nElementsOnThisProc();
//  Epetra_Map Map(-1, locElements, 0, *m_commPtr);

  bool isSolverInterfaceSet = m_solverInterface->isSetupForStabilityRun();
  if (isSolverInterfaceSet != true)
  {
    if (m_commPtr->MyPID () == 0)
    {
      cout << "StabilityEveluator::computeDominantModes returned with error. solverInterface is not setup for stability run" << endl;
    }
    return 1;
  }

  Epetra_Map Map = m_solverInterface->getEpetraMap(m_commPtr);

  RCP<MV> ivec = rcp (new MV (Map, a_blockSize, m_solverInterface) );
  MVT::MvRandom( *ivec );

  RCP<OP> Aop = rcp (new OP ( m_solverInterface ) );

  RCP<Anasazi::BasicEigenproblem<double,BaseMV,BaseOP> > EigenProblem = rcp ( new Anasazi::BasicEigenproblem<double,BaseMV,BaseOP>(Aop, ivec) );

  EigenProblem->setHermitian(a_isOpSymmetric);
  EigenProblem->setNEV( a_nev );
    
  bool boolret = EigenProblem->setProblem();
  if (boolret != true)
  {
    if (m_commPtr->MyPID () == 0)
    {
      cout << "Anasazi::BasicEigenproblem::setProblem() returned with error." << endl;
    }
    return 1;
  }
  else
  {
    if (m_commPtr->MyPID () == 0)
    {
       cout << "Anasazi::BasicEigenproblem::setProblem() returned with success." << endl;
    }
  }

  // Initialize Block Arnoldi solver
  Anasazi::BlockKrylovSchurSolMgr<double, BaseMV, BaseOP> EigenSolverMgr(EigenProblem, EigenPL);

  // Solve the problem to the specied tolerances or length
  Anasazi::ReturnType returnCode = EigenSolverMgr.solve();
  if (returnCode != Anasazi::Converged)
  {
    if (m_commPtr->MyPID() == 0)
    {
      cout << "Anasazi::EigensolverMgr::solve() returned unconverged." << endl;
    }
    return 1;
  }
  else
  {
    if (m_commPtr->MyPID() == 0)
    {
      cout << "Anasazi::EigensolverMgr::solve() returned converged." << endl;
    }
  }

  // Get Ritz values from EigenSolver
  std::vector<Anasazi::Value<double> > ritzValues = EigenSolverMgr.getRitzValues();

  // Output computed eigenvalues and their direct residuals
  if (a_verbose && m_commPtr->MyPID()==0)
  {

    int numritz = (int)ritzValues.size();
    std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
    std::cout<<std::endl<< "Computed Ritz Values"<< std::endl;
    if (EigenProblem->isHermitian()) {
      std::cout<< std::setw(16) << "Real Part"
        << std::endl;
      std::cout<<"-----------------------------------------------------------"<<std::endl;
      for (int i=0; i<numritz; i++) {
        std::cout<< std::setw(16) << ritzValues[i].realpart
          << std::endl;
      }
      std::cout<<"-----------------------------------------------------------"<<std::endl;
    }
    else {
      std::cout<< std::setw(16) << "Real Part"
        << std::setw(16) << "Imag Part"
        << std::endl;
      std::cout<<"-----------------------------------------------------------"<<std::endl;
      for (int i=0; i<numritz; i++) 
      {
        std::cout<< std::setw(16) << ritzValues[i].realpart
          << std::setw(16) << ritzValues[i].imagpart
          << std::endl;
      }
      std::cout<<"-----------------------------------------------------------"<<std::endl;
    }
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  Anasazi::Eigensolution<double,BaseMV> sol = EigenProblem->getSolution ();
  std::vector<Anasazi::Value<double> > evals = sol.Evals;
  RCP<BaseMV> evecs = sol.Evecs;
  std::vector<int> index = sol.index;
  int numev = sol.numVecs;

  if (numev > 0)
  {
    // plot eigenvectors
    if (a_plotEVComps.size() > 0)
    {
      Anasazi::EpetraMultiVecSolverExt* castEigVec = dynamic_cast<Anasazi::EpetraMultiVecSolverExt*>(evecs.get()); 
      TEUCHOS_TEST_FOR_EXCEPTION( castEigVec==NULL,  std::invalid_argument, "StabilityEvaluator::computeDominantModes cast of MultiVec<double> to EpetraMultiVecSolverExt failed.");

      for (int icomp = 0; icomp < a_plotEVComps.size(); icomp++)
      {
        std::string fileName = "computed_evec_comp_" + SSTR(icomp);
        m_solverInterface->plotEpetraVector(*((*castEigVec)(icomp)), fileName);
      }
    }

    // Compute residuals
    Teuchos::LAPACK<int,double> lapack;
    std::vector<double> normA(numev);

    if (EigenProblem->isHermitian())
    {
      // Get storage
      MV Aevecs(Map, numev, m_solverInterface);
      Teuchos::SerialDenseMatrix<int,double> B(numev, numev);
      B.putScalar(0.0);
      for (int i=0; i<numev; i++)
      {
        B(i,i) = evals[i].realpart;
      }

      // Compute A*evecs
      OPT::Apply( *Aop, *evecs, Aevecs );

      // Compute A*evecs - lambda*evecs and its norm
      MVT::MvTimesMatAddMv( -1.0, *evecs, B, 1.0, Aevecs);
      MVT::MvNorm( Aevecs, normA );

      // Scale the norms by the eigenvalue
      for (int i = 0; i < numev; i++)
      {
        normA[i] /= Teuchos::ScalarTraits<double>::magnitude( evals[i].realpart);
      }
    }
    else
    {
      // The problem in non-Hermitian.
      int i = 0;
      std::vector<int> curind(1);
      std::vector<double> resnorm(1), tempnrm(1);
      RCP<BaseMV> tempAevec;
      RCP<const BaseMV> evecr, eveci;
      MV Aevec(Map, numev, m_solverInterface);

      // Compute A*evecs
      OPT::Apply( *Aop, *evecs, Aevec );

      Teuchos::SerialDenseMatrix<int,double> Breal(1,1), Bimag(1,1);
      while (i < numev)
      {  
        if (index[i] == 0)
        {
          // Get a view of the current eigenvector (evecr)
          curind[0] = i;
          evecr = MVT::CloneView( *evecs, curind );

          // Get a copy of A*evecr
          tempAevec = MVT::CloneCopy( Aevec, curind );

          // Compute A*evecr - lambda*evecr
          Breal(0,0) = evals[i].realpart;
          MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );

          // Compute the norm of the residual and increment counter
          MVT::MvNorm( *tempAevec, resnorm );
          normA[i] = resnorm[0]/Teuchos::ScalarTraits<MagnitudeType>::magnitude( evals[i].realpart );
          i++;
        }  
        else
        {
          // Get a view of the real part of the eigenvector (evecr)
          curind[0] = i;
          evecr = MVT::CloneView( *evecs, curind );
          // Get a copy of A*evecr
          tempAevec = MVT::CloneCopy( Aevec, curind );
          // Get a view of the imaginary part of the eigenvector (eveci)
          curind[0] = i+1;
          eveci = MVT::CloneView( *evecs, curind );
          // Set the eigenvalue into Breal and Bimag
          Breal(0,0) = evals[i].realpart;
          Bimag(0,0) = evals[i].imagpart;
          // Compute A*evecr - evecr*lambdar + eveci*lambdai
          MVT::MvTimesMatAddMv( -1.0, *evecr, Breal, 1.0, *tempAevec );
          MVT::MvTimesMatAddMv( 1.0, *eveci, Bimag, 1.0, *tempAevec );
          MVT::MvNorm( *tempAevec, tempnrm );
          // Get a copy of A*eveci
          tempAevec = MVT::CloneCopy( Aevec, curind );
          // Compute A*eveci - eveci*lambdar - evecr*lambdai
          MVT::MvTimesMatAddMv( -1.0, *evecr, Bimag, 1.0, *tempAevec );
          MVT::MvTimesMatAddMv( -1.0, *eveci, Breal, 1.0, *tempAevec );
          MVT::MvNorm( *tempAevec, resnorm );
          // Compute the norms and scale by magnitude of eigenvalue
          normA[i] = lapack.LAPY2( tempnrm[0], resnorm[0] ) /
            lapack.LAPY2( evals[i].realpart, evals[i].imagpart );
          normA[i+1] = normA[i];
          i=i+2;
        }
      }
    }

    // Output computed eigenvalues and their direct residuals
    if (a_verbose && m_commPtr->MyPID()==0) 
    {
      std::cout.setf(std::ios_base::right, std::ios_base::adjustfield);
      std::cout<<std::endl<< "Actual Residuals"<<std::endl;
      if (EigenProblem->isHermitian()) 
      {
        std::cout<< std::setw(16) << "Real Part"
          << std::setw(20) << "Direct Residual"<< std::endl;
        std::cout<<"-----------------------------------------------------------"<<std::endl;
        for (int i=0; i<numev; i++) 
        {
          std::cout<< std::setw(16) << evals[i].realpart
            << std::setw(20) << normA[i] << std::endl;
        }
        std::cout<<"-----------------------------------------------------------"<<std::endl;
      }
      else 
      {
        std::cout<< std::setw(16) << "Real Part"
          << std::setw(16) << "Imag Part"
          << std::setw(20) << "Direct Residual"<< std::endl;
        std::cout<<"-----------------------------------------------------------"<<std::endl;
        for (int i=0; i<numev; i++) 
        {
          std::cout<< std::setw(16) << evals[i].realpart
            << std::setw(16) << evals[i].imagpart
            << std::setw(20) << normA[i] << std::endl;
        }
        std::cout<<"-----------------------------------------------------------"<<std::endl;
      }
    }
  } 

  return 0;
}

/*********/
