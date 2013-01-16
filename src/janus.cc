/*
Copyright (c) 2012, Bruno Turcksin.

This file is part of Janus.

Janus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
he Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Janus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Janus.  If not, see <http://www.gnu.org/licenses/>.
*/

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
#include "DIFFUSION_SOLVER.hh"
#include "EXCEPTION.hh"
#include "TRANSPORT_SOLVER.hh"

using namespace std;

int main(int argc,char** argv)
{
  try
  {
    Check(argc==6,"Wrong number of inputs");

    string geometry_filename(argv[2]);
    string parameters_filename(argv[3]);
    string cross_sections_filename(argv[4]);
    string output_filename(argv[5]);

    MPI_Init(&argc,&argv);
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    if (strcmp(argv[1],"--t")==0)
    {
      cout<<"Initialization"<<endl;
      TRANSPORT_SOLVER transport_solver(&geometry_filename,&parameters_filename,
          &cross_sections_filename,&output_filename,&comm);
      cout<<"Start solving"<<endl;
      transport_solver.Solve();
      cout<<"End solving"<<endl;
      cout<<"Start writing output file"<<endl;
      transport_solver.Write_in_file();
      cout<<"End writing output file"<<endl;
    }
    else
    {
      Check(strcmp(argv[1],"--d")==0,string ("Unknown option: only options available are --t and --d."));
      cout<<"Initialization"<<endl;
      DIFFUSION_SOLVER diffusion_solver(&geometry_filename,&parameters_filename,
          &cross_sections_filename,&output_filename,&comm);
      cout<<"Start solving"<<endl;
      diffusion_solver.Solve();
      cout<<"End solving"<<endl;
      cout<<"Start writing output file"<<endl;
      diffusion_solver.Write_in_file();
      cout<<"End writing output file"<<endl;
    }
  }
  catch (const string error)
  {
    cerr<<error<<endl;
    exit(1);
  }

  MPI_Finalize();

  return 0;
}
