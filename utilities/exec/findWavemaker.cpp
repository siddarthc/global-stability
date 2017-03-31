/*
*
*
*
*
*/

#include "ParmParse.H"
#include "parstream.H"

#include "EBIndexSpace.H"
#include "EBISLayout.H"
#include "BoxIterator.H"
#include "RealVect.H"

#include "AMRINSUtils.H"

#include <iostream>
#include "memusage.H"
#include "memtrack.H"

#include "WavemakerUtils.H"

#include "UsingNamespace.H"

/********/
int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif

  { //scoping trick

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

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

    ProblemDomain coarsestDomain;
    Real coarsestDx;

    // specializing for INS:
    AMRParameters params;
    getAMRINSParameters(params, coarsestDomain);

    //define geometry
    AMRINSGeometry(params, coarsestDomain);

    int nModes;
    pp.get("num_modes", nModes);

    Vector<std::string> fileName;
    pp.getarr("file_name", fileName, 0, nModes);
    Vector<std::string> directModeFile;
    pp.getarr("direct_mode_file", directModeFile, 0, nModes);
    Vector<std::string> adjointModeFile;
    pp.getarr("adjoint_mode_file", adjointModeFile, 0, nModes);

    coarsestDx = params.m_domainLength/Real(coarsestDomain.size(0));

    for (int iMode = 0; iMode < nModes; iMode++)
    {
      WavemakerUtils::plotWavemaker(fileName[iMode], directModeFile[iMode], adjointModeFile[iMode], params.m_refRatio, coarsestDomain, coarsestDx);
    }

  } // end scoping

#ifdef CH_MPI
  CH_TIMER_REPORT();
  dumpmemoryatexit();
  MPI_Finalize();
#endif
}
