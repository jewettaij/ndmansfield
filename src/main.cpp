// ---- Compiling instructions: ----
//
//
//
// To compile for profiling, use:
// rm gmon.out a.out ndmansfield ndmansfield.o
// g++ -c -O4 -DNDEBUG -funroll-loops -ffast-math ndmansfield.cc -pg
// g++ -O4 -DNDEBUG -funroll-loops -ffast-math -o ndmansfield ndmansfield.o -pg
// gprof ndmansfield > gprof_results.txt
// For more information on profiling:
// http://www.cs.utah.edu/dept/old/texinfo/as/gprof.html


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>     // defines stringstream
#include <cstdlib>     // defines atod()
#include <cstring>     // defines strcmp()
using namespace std;

#include "ndmansfield.h"   // Contains a description of the NDmansfield data structure
                        // we will be using to store the state of our system.

#include "random_gen.h" // Some random number generators.  This file defines
                        // the functions RANDOM_REAL_0_1() and INIT_RANDOM()

#include "hamiltonian.h" // The object that will calculate the energy.

#include "analysis.h"   // Diagnostics to report while running the simulation.




void ReadFile(string fname, NDmansfield& ndmansfield); // defined below






string g_program_name("ndmansfield");
string g_version_string("0.11");
string g_date_string("2016-3-29");






int main(int argc, char **argv)
{
  stringstream err_msg;
  string explanation_of_syntax("Usage Syntax: \n"
                               "ndmansfield -box xsize ysize zsize -tsave tsave [options] > trajectory.raw");

  cerr << g_program_name << " v" << g_version_string << ", " 
       << g_date_string << "\n" << flush;

  g_dim = -1;
  const long long DISABLE = -1;
  long long tstart = 0;
  long long tstop = DISABLE;
  long long tsave = DISABLE;
  long seed = -1;
  long *aSize = NULL;
  string starting_filename;
  bool record_cyclic = false;
  double bend_energy_coeff = 0.0;
  double twist_energy_coeff = 0.0;
  
  try {
    //Read in the command line argument(s) passed by the user.

    bool syntax_error_occured = false;
    long i=1;
    while (i < argc)
    {
      int ndelete = 1;
      if ((strncmp(argv[i],"-h",2)==0) ||
          (strncmp(argv[i],"--h",3)==0) ||
          (strncmp(argv[i],"--help",3)==0) ||
          (strncmp(argv[i],"-H",2)==0) ||
          (strncmp(argv[i],"--H",3)==0) ||
          (strncmp(argv[i],"-?",2)==0) ||
          (strncmp(argv[i],"--?",3)==0))
      {
        cerr << explanation_of_syntax;
        exit(0);
      }
      else if (strcmp(argv[i], "-tsave")==0) {
        if ((i+1>=argc) || (strlen(argv[i+1])==0) || (! isdigit(argv[i+1][0])))
        {
          err_msg << "Error: Expected an integer following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          tsave = atoll(argv[i+1]);
        ndelete = 2;
      }
      else if (strcmp(argv[i], "-tstart")==0) {
        if ((i+1>=argc) || (strlen(argv[i+1])==0) || (! isdigit(argv[i+1][0])))
        {
          err_msg << "Error: Expected an integer following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          tstart = atoll(argv[i+1]);
        ndelete = 2;
      }
      else if (strcmp(argv[i], "-tstop")==0) {
        if ((i+1>=argc) || (strlen(argv[i+1])==0) || (! isdigit(argv[i+1][0])))
        {
          err_msg << "Error: Expected a number following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          tstop = atoll(argv[i+1]);
        ndelete = 2;
      }
      else if (strcmp(argv[i], "-seed")==0) {
        if ((i+1>=argc) || (strlen(argv[i+1])==0) || (! isdigit(argv[i+1][0])))
        {
          err_msg << "Error: Expected a number following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          seed = atoll(argv[i+1]);
        ndelete = 2;
      }
      else if (strcmp(argv[i], "-startcrd")==0) {
        if ((i+1 >= argc) || (strlen(argv[i+1]) == 0))
        {
          err_msg << "Error: Expected a file name following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          starting_filename.assign(argv[i+1]);
        ndelete = 2;
      }
      else if (strcmp(argv[i], "-cyclic")==0) {
        if ((i+1 >= argc) || (strlen(argv[i+1]) == 0))
        {
          err_msg << "Error: Expected either \"yes\" or \"no\" following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          record_cyclic = (strcmp(argv[i+1], "yes") == 0);
        ndelete = 2;
      }
      else if (strcmp(argv[i], "-box")==0)
      {
        //------------------
        //If "g_dim" is a constant, then comment out the next 11 lines:

        g_dim = 0;
        while ((i+g_dim+1 < argc) && 
               (strlen(argv[i+g_dim+1])>0) && 
               isdigit(argv[i+g_dim+1][0])) {
          g_dim++;
        }

        if (g_dim < 2) {
          err_msg << "Error: Expected at least 2 integers following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        aSize = new long [g_dim];
        for (int d=0; d < g_dim; d++) {
          //------------------
          // ...and uncomment the next 7 lines:
          //if ((strlen(argv[i+1]) == 0) || (! isdigit(argv[i+1][0]))) {
          //    err_msg << "Error: Expected " << g_dim 
          //            << "\" integers following the \""
          //            << argv[i] << "\" argument.\n";
          //    throw InputErr(err_msg.str());
          //}
          //else
          aSize[d] = atoi(argv[i+1+d])+2;
        }
        ndelete = 1 + g_dim;
      }
      else if ((strcmp(argv[i], "-bend-energy")==0) ||
               (strcmp(argv[i], "-ebend")==0))
      {
        if ((i+1>=argc) || (strlen(argv[i+1])==0))
        {
          err_msg << "Error: Expected a number following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          bend_energy_coeff = atof(argv[i+1]);
        ndelete = 2;
      }
      else if ((strcmp(argv[i], "-twist-energy")==0) ||
               (strcmp(argv[i], "-etwist")==0) ||
               (strcmp(argv[i], "-dihedral-energy")==0) ||
               (strcmp(argv[i], "-edihdral")==0) ||
               (strcmp(argv[i], "-edih")==0))
      {
        if ((i+1>=argc) || (strlen(argv[i+1])==0))
        {
          err_msg << "Error: Expected a number following the \""
                  << argv[i] << "\" argument.\n";
          throw InputErr(err_msg.str());
        }
        else
          twist_energy_coeff = atof(argv[i+1]);
        ndelete = 2;
      }
      else {
          err_msg << "Error: Unrecognized argument: \""
                  << argv[i] << "\"\n";
          throw InputErr(err_msg.str());
      }

      i += ndelete;

    } // while (i < argc)
    if (g_dim < 2) {
      string syntax_example =
        "Usage Syntax: \n"
        "\n"
        "   ndmansfield -box xsize ysize zsize -tsave tsave [options] > trajectory.raw\n"
        "\n"
        "   (Note: The number of arguments following the -box command varies\n"
        "    depending on the number of dimensions your polymer lives in.\n"
        "    In 2-D, you would omit the \"zsize\" argument, for example.\n"
        "    See the \"docs_ndmansfield.txt\" file for detailed instructions.)\n";
      err_msg << "Error: Missing \"-box\" command.\n\n"
              << syntax_example << endl;
      throw InputErr(err_msg.str());
    }
    assert(aSize);
  }
  catch (InputErr& e) {
    cerr << e.what() << endl;
    exit(ERR_SYNTAX);
  }



  try {

    if (tsave != DISABLE)
      cerr << " backup interval: " << tsave << "\n"
        "        (Note: Coordinates will be printed to the terminal (stdout)\n"
        "               Use \"> file_name.raw\" to redirect this to a file.)\n";


    else
      cerr <<
        " ####################################################################\n"        " #######      WARNING: NO COORDINATE DATA WILL BE SAVED.      #######\n"
        " #######         ARE YOU SURE THIS IS WHAT YOU WANT?          #######\n"
        " ####### (Use the \"-tsave\" argument to specify how frequently #######\n"
        " #######     to print coordinates to the terminal/stdout.)    #######\n"
        " ####################################################################\n";
        
    cerr << " record_only_cyclic: " << (record_cyclic ? "true\n" : "false\n");
    if (tstop != DISABLE)
      cerr << " simulation tstop: "<< tstop <<" MC moves\n";
    cerr << "   (starting time: t=" << tstart << ")\n"
         << " number of dimensions: " << g_dim << "\n"
         << " lattice size: ";
    //for(int d=0; d<g_dim; ++d) {
    //  cerr << aSize[d];
    //  if (d < g_dim-1) cerr << "x";
    //}
    //cerr << "\n"
    //     << " available lattice size: ";
    for(int d=0; d<g_dim; ++d) {
      cerr << aSize[d]-2;
      if (d < g_dim-1) cerr << "x";
    }
    cerr << "\n";
    if (starting_filename != "")
      cerr << " Initial coordinates in: \"" << starting_filename << "\"\n" << flush;




    //Setup the random-number generator:
    if (seed == -1)
      RANDOM_INIT();
    else
      RANDOM_INIT(seed);




    //Create the objects we will need:

    NDmansfield ndmansfield(aSize);
    if (starting_filename != "")
      ReadFile(starting_filename, ndmansfield);

    Hamiltonian hamiltonian(ndmansfield, bend_energy_coeff, twist_energy_coeff);

    double energy = hamiltonian.CalcEnergy(ndmansfield);

    long long num_moves_accepted = 0; //a "crazy" initial value




    // ---- Run the simulation: -----

    // Main loop:

    long long t_prev = 0;
    for(long long t=tstart; 
        ((t <= tstop) || (tstop == DISABLE)); 
        ++t)
    {
      bool record_this = true;
      if (tsave == DISABLE)
        record_this = false;
      if ((t-t_prev) < tsave)
        record_this = false;
      if ((record_cyclic) && (! ndmansfield.IsCyclic()))
        record_this = false;
      if ((t == tstart) && (starting_filename == ""))
        record_this = true;
      if (record_this) {

        cerr << "   t = " << t << endl;

        // Calculate various metrics of the chain and report them to the user.
        // This should help them determine whether then chain has equilibrated,
        // (and whether or not it is too soon to halt the simulation).
        double aSpatialCorrelation[g_dim];
        SpatialCorrelation(aSpatialCorrelation,
                           ndmansfield);

        cerr << "      correlation_of_last_coordinate_vs_i = ";
        //for (int d=0; d+1<g_dim; d++)
        //  cerr << aSpatialCorrelation[d] << " ";
        cerr << aSpatialCorrelation[g_dim-1] << endl;

        long aBondCountPlusXYZ[g_dim];
        long aBondCountMinusXYZ[g_dim];
        CountBondDirections(aBondCountPlusXYZ,
                            aBondCountMinusXYZ,
                            ndmansfield);
        
        cerr << "      bond_count_xyz = ";
        for (int d=0; d+1<g_dim; d++)
          cerr << (aBondCountPlusXYZ[d] + aBondCountMinusXYZ[d]) << " ";
        cerr << (aBondCountPlusXYZ[g_dim-1] + aBondCountMinusXYZ[g_dim-1]) << endl;

        cout << ndmansfield; // write the state of the system to a file.
        cout << flush;
        t_prev = t;
      }
      //cerr.precision(11);
      //cerr << t << " " << energy << endl;

      // Advance the dynamics of the system...

      if (ndmansfield.MonteCarloStep(hamiltonian, energy))
        num_moves_accepted++;
    } // for(long long t=0; t < tstop; ++t)

    cerr << "   Evolved until t = " << tstop << "\n"
         << "  exiting...\n" << flush;

  }
  catch (InputErr& e) {
    cerr << e.what() << endl;
    exit(ERR_INPUT);
  }


} // main()





void ReadFile(string fname, NDmansfield& ndmansfield)
{
  try {    
    ifstream
      trajectory_file(fname.c_str(), ios::in);
    if (! trajectory_file) {
      stringstream err_msg;
      err_msg << "Error: unable to open file \"" 
              << fname << "\" for reading" << endl;
      throw InputErr(err_msg.str());
    }
    trajectory_file >> ndmansfield;
  }
  catch (InputErr& e) {
    cerr << e.what() << endl;
    exit(ERR_INPUT);
  }
} //ReadFile()


