#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include "mpi.h"
#include "gsl_math.h"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_MultiVector.h"
#include "CELL.hh"
#include "EDGE.hh"
#include "LS.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
#include "TRANSPORT_OPERATOR.hh"
#include "TRIANGULATION.hh"

using namespace std;

typedef vector<double> d_vector;

int main(int argc,char** argv)
{
  unsigned int n_sources(1);
  unsigned int n_materials(1);
  unsigned int flux_moments_size(0);
  const double sqrt_4pi(2.*sqrt(M_PI));
  string geometry_inp("/home/bruno/Documents/Transport/janus/tests/geometry_transport_1.inp");
  string parameters_inp("/home/bruno/Documents/Transport/janus/tests/parameters_transport_1.inp");

  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // Reference solution
  d_vector solution(36,0.);
  solution[0] = sqrt_4pi * 0.974449442413;
  solution[1] = sqrt_4pi * 0.97444944203;
  solution[2] = sqrt_4pi * 1.66326834888;
  solution[3] = sqrt_4pi * 1.66326833258;
  solution[4] = sqrt_4pi * 0.974449442389;
  solution[5] = sqrt_4pi * 0.974449442389;
  solution[6] = sqrt_4pi * 1.66326834662;
  solution[7] = sqrt_4pi * 1.66326834662;
  solution[8] = sqrt_4pi * 0.97444944209;
  solution[9] = sqrt_4pi * 0.974449442413;
  solution[10] = sqrt_4pi * 1.66326833258;
  solution[11] = sqrt_4pi * 1.66326834888;
  solution[12] = sqrt_4pi * 1.66833181666;
  solution[13] = sqrt_4pi * 1.6683318057;
  solution[14] = sqrt_4pi * 1.6683318057;
  solution[15] = sqrt_4pi * 1.66833181666;
  solution[16] = sqrt_4pi * 1.66833180693;
  solution[17] = sqrt_4pi * 1.66833180693;
  solution[18] = sqrt_4pi * 1.66833180693;
  solution[19] = sqrt_4pi * 1.66833180693;
  solution[20] = sqrt_4pi * 1.6683318057;
  solution[21] = sqrt_4pi * 1.66833181666;
  solution[22] = sqrt_4pi * 1.66833181666;
  solution[23] = sqrt_4pi * 1.6683318057;
  solution[24] = sqrt_4pi * 1.66326833258;
  solution[25] = sqrt_4pi * 1.66326834888;
  solution[26] = sqrt_4pi * 0.97444944209;
  solution[27] = sqrt_4pi * 0.974449442413;
  solution[28] = sqrt_4pi * 1.66326834662;
  solution[29] = sqrt_4pi * 1.66326834662;
  solution[30] = sqrt_4pi * 0.974449442389;
  solution[31] = sqrt_4pi * 0.974449442389;
  solution[32] = sqrt_4pi * 1.66326834888;
  solution[33] = sqrt_4pi * 1.66326834888;
  solution[34] = sqrt_4pi * 0.974449442413;
  solution[35] = sqrt_4pi * 0.97444944209;

  TRIANGULATION triangulation(&geometry_inp);
  PARAMETERS parameters(&parameters_inp);
  vector<QUADRATURE*> quad(1,NULL);

  // Create the triangulation
  triangulation.Read_geometry();
  triangulation.Build_edges();

  // Create the parameters
  parameters.Read_parameters(n_sources,n_materials);

  // Build the quadrature
  quad[0] = new LS(parameters.Get_sn_order(),parameters.Get_L_max(),false);
  quad[0]->Build_quadrature();

  // Build the dof handler
  DOF_HANDLER dof_handler(&triangulation,parameters);
  dof_handler.Compute_sweep_ordering(quad);

  // Flux moments map and vector
  flux_moments_size = dof_handler.Get_n_dof()*quad[0]->Get_n_mom() +
    dof_handler.Get_n_sf_per_dir()*quad[0]->Get_n_dir();
  Epetra_Map flux_moments_map(flux_moments_size,0,comm);
  Epetra_MultiVector flux_moments(flux_moments_map,1);

  // Solve the transport
  TRANSPORT_OPERATOR transport_operator(&dof_handler,&parameters,quad[0],&comm,
      &flux_moments_map);

  // Compute right-hand side for BiCGSTAB
  Epetra_MultiVector rhs(flux_moments);
  transport_operator.Sweep(rhs,true);

  Epetra_LinearProblem problem(&transport_operator,&flux_moments,&rhs);

  AztecOO solver(problem);

  // Set Aztec parameters
  solver.SetAztecOption(AZ_precond,AZ_none);
  solver.SetAztecOption(AZ_solver,AZ_bicgstab);
  solver.SetAztecOption(AZ_conv,AZ_rhs);

  // Solve the transport equation
  solver.Iterate(parameters.Get_max_it(),parameters.Get_tolerance());

  // Apply the preconditioner to get the solution
  MIP* precond(transport_operator.Get_mip());
  precond->Solve(flux_moments);

  for (unsigned int i=0; i<dof_handler.Get_n_dof(); ++i)
    assert(fabs(flux_moments[0][i]-solution[i])<0.0001);
  
  delete quad[0];
  quad[0] = NULL;

  MPI_Finalize();

  return 0;
}