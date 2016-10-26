/*
 *
 *
 *
 *
 */

#include "ParmParse.H"
#include "parstream.H"
#include "memusage.H"
#include "memtrack.H"
#include "CH_Attach.H"
#include "CH_HDF5.H"
#include "StabilityEvaluator.H"

// Trilinos library includes
#ifdef CH_MPI
#include <Epetra_MpiComm.h>
#else
#include <Eperta_SerialComm.h>
#endif

/*********/
void
exampleRoutine (const Epetra_Comm& comm,
                std::ostream& out)
{
  if (comm.MyPID () == 0) {
    // On (MPI) Process 0, print out the Epetra software version.
    out << "in main:" << std::endl << std::endl;
  }
}

/*********/
int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif

  { //scoping trick

 // Check for an input file
  char* inFile = NULL;

  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      pout() << "Usage: <executable name> <inputfile>" << endl;
      pout() << "No input file specified" << endl;
      return -1;
    }
  
  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

  }

#ifdef CH_MPI

  // Trilinos MPI test
  MPI_Comm wcomm = Chombo_MPI::comm;
/*
  Epetra_MpiComm comm (wcomm);

  const int myRank = comm.MyPID();
  const int numProcs = comm.NumProc();

  if (myRank == 0)
    {
      std::cout << "Total number of processes: " << numProcs << std::endl;
    }

  exampleRoutine(comm, std::cout);
*/
  StabilityEvaluator<double> testEvaluator(1,1.0,1.0,"dummy",&wcomm);
  testEvaluator.exampleRoutine();
/*
  if (myRank == 0)
    {
      std::cout << "End Result: TEST PASSED" << std::endl;
    } 
*/
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

}
