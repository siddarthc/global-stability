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

  StabilityEvaluator<Real> evaluator;

  }

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif

}
