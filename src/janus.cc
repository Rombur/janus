/**
 * \mainpage
 * Documentation of Janus. This code solves the one-group transport equation
 * for neutrons on rectangular domains using the Sn method and discontinuous
 * finite elements. The input files of the code can be generated using Diana
 * (geometry) and Mercury (parameters). The output file can be converted to a
 * silo file using Apollo.
 */


#include <string>
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "EXCEPTION.hh"
#include "TRANSPORT_SOLVER.hh"

using namespace std;

int main(int argc,char** argv)
{
  Check(argc==4,"Wrong number of inputs");

  string geometry_filename(argv[1]);
  string parameters_filename(argv[2]);
  string output_filename(argv[3]);

  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  try
  {
    cout<<"Initialization"<<endl;
    TRANSPORT_SOLVER transport_solver(&geometry_filename,&parameters_filename,
        &output_filename,&comm);
    cout<<"Start solving"<<endl;
    transport_solver.Solve();
    cout<<"End solving"<<endl;
    cout<<"Start writing output file"<<endl;
    transport_solver.Write_in_file();
    cout<<"End writing output file"<<endl;
  }
  catch (const string error)
  {
    cerr<<error<<endl;
    exit(1);
  }

  MPI_Finalize();

  return 0;
}
