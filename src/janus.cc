#include <cassert>
#include <string>
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "TRANSPORT_SOLVER.hh"

using namespace std;

int main(int argc,char** argv)
{
  assert(argc==4);

  string geometry_filename(argv[1]);
  string parameters_filename(argv[2]);
  string output_filename(argv[3]);

  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  cout<<"Initialization"<<endl;
  TRANSPORT_SOLVER transport_solver(&geometry_filename,&parameters_filename,
      &output_filename,&comm);
  cout<<"Start solving"<<endl;
  transport_solver.Solve();
  cout<<"End solving"<<endl;
  cout<<"Start writing output file"<<endl;
  transport_solver.Write_in_file();
  cout<<"End writing output file"<<endl;

  MPI_Finalize();

  return 0;
}
