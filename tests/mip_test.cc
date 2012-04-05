#include <cassert>
#include <cmath>
#include <string>
#include "gsl_math.h"
#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "DOF_HANDLER.hh"
#include "GLC.hh"
#include "MIP.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
#include "TRIANGULATION.hh"

using namespace std;

int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const double inv_sqrt_4pi(1./sqrt((2.*M_PI)));

  string geometry_inp("/home/bruno/Documents/Transport/janus/tests/geometry_mip.inp");
  string parameters_inp("/home/bruno/Documents/Transport/janus/tests/parameters_mip.inp");

  TRIANGULATION triangulation(&geometry_inp);
  PARAMETERS parameters(&parameters_inp);
  vector<QUADRATURE*> quad(1,NULL);
   
  triangulation.Read_geometry();
  triangulation.Build_edges();

  parameters.Read_parameters(triangulation.Get_n_sources(),
      triangulation.Get_n_materials());
  
  quad[0] = new GLC(parameters.Get_sn_order(),parameters.Get_L_max(),
      parameters.Get_galerkin());
  quad[0]->Build_quadrature();

  DOF_HANDLER dof_handler(&triangulation,parameters);
  dof_handler.Compute_sweep_ordering(quad);

  Epetra_Map map(dof_handler.Get_n_dof(),0,comm);
  Epetra_MultiVector flux_moments(map,1);
  Epetra_MultiVector rhs(map,1);
  Epetra_MultiVector solution(map,1);
  
  flux_moments[0][0] = 0.15480683;
  flux_moments[0][1] = 0.22424409;
  flux_moments[0][2] = 0.354209;
  flux_moments[0][3] = 0.22424409;
  flux_moments[0][4] = 0.22424409;
  flux_moments[0][5] = 0.15480683;
  flux_moments[0][6] = 0.22424409;
  flux_moments[0][7] = 0.354209;
  flux_moments[0][8] = 0.22424409;
  flux_moments[0][9] = 0.354209;
  flux_moments[0][10] = 0.22424409;
  flux_moments[0][11] = 0.15480683;
  flux_moments[0][12] = 0.354209;
  flux_moments[0][13] = 0.22424409;
  flux_moments[0][14] = 0.15480683;
  flux_moments[0][15] = 0.22424409;

  solution[0][0] = inv_sqrt_4pi * 0.01579434;
  solution[0][1] = inv_sqrt_4pi * 0.0410164;
  solution[0][2] = inv_sqrt_4pi * 0.11759913;
  solution[0][3] = inv_sqrt_4pi * 0.04101419;
  solution[0][4] = inv_sqrt_4pi * 0.04100877;
  solution[0][5] = inv_sqrt_4pi * 0.01578621;
  solution[0][6] = inv_sqrt_4pi * 0.04109177;
  solution[0][7] = inv_sqrt_4pi * 0.11758109;
  solution[0][8] = inv_sqrt_4pi * 0.04103222;
  solution[0][9] = inv_sqrt_4pi * 0.11757449;
  solution[0][10] = inv_sqrt_4pi * 0.04109206;
  solution[0][11] = inv_sqrt_4pi * 0.01576069;
  solution[0][12] = inv_sqrt_4pi * 0.1176016;
  solution[0][13] = inv_sqrt_4pi * 0.04107562;
  solution[0][14] = inv_sqrt_4pi * 0.01581032;
  solution[0][15] = inv_sqrt_4pi * 0.04107403;

  MIP precond(&dof_handler,&parameters,quad[0],&comm);
  precond.Solve(flux_moments);

  for (unsigned int i=0; i<dof_handler.Get_n_dof(); ++i)
    assert(fabs(flux_moments[0][i]-solution[0][i])<0.001);

  delete quad[0];
  quad[0] = NULL;

  MPI_Finalize();

  return 0;
}
