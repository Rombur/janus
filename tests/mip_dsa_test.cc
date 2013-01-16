#include <cassert>
#include <cmath>
#include <string>
#include "gsl_math.h"
#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "CROSS_SECTIONS.hh"
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

  const double four_pi(4.*M_PI);
  string cross_sections_inp("cross_sections_mip_dsa.inp");
  string geometry_inp("geometry_mip_dsa.inp");
  string parameters_inp("parameters_mip_dsa.inp");

  CROSS_SECTIONS cross_sections(&cross_sections_inp);
  TRIANGULATION triangulation(&geometry_inp);
  PARAMETERS parameters(&parameters_inp);
  vector<QUADRATURE*> quad(1,NULL);
   
  triangulation.Read_geometry();
  triangulation.Build_edges();

  parameters.Read_parameters(triangulation.Get_n_sources());

  cross_sections.Read_regular_cross_sections(triangulation.Get_n_materials(),
      parameters.Get_permutation_type(),false);
  cross_sections.Apply_ang_lvls_and_tc(parameters.Get_multigrid(),
      parameters.Get_transport_correction(),parameters.Get_optimal_tc(),
      triangulation.Get_n_materials(),parameters.Get_sn_order());
  
  quad[0] = new GLC(parameters.Get_sn_order(),cross_sections.Get_L_max(),
      parameters.Get_galerkin());
  quad[0]->Build_quadrature(four_pi);

  DOF_HANDLER dof_handler(&triangulation,parameters,cross_sections);
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

  solution[0][0] = 0.01579434;
  solution[0][1] = 0.0410164;
  solution[0][2] = 0.11759913;
  solution[0][3] = 0.04101419;
  solution[0][4] = 0.04100877;
  solution[0][5] = 0.01578621;
  solution[0][6] = 0.04109177;
  solution[0][7] = 0.11758109;
  solution[0][8] = 0.04103222;
  solution[0][9] = 0.11757449;
  solution[0][10] = 0.04109206;
  solution[0][11] = 0.01576069;
  solution[0][12] = 0.1176016;
  solution[0][13] = 0.04107562;
  solution[0][14] = 0.01581032;
  solution[0][15] = 0.04107403;

  MIP precond(&comm,&parameters,&dof_handler,quad[0]);
  precond.Solve(flux_moments);

  for (unsigned int i=0; i<dof_handler.Get_n_dof(); ++i)
    assert(fabs(flux_moments[0][i]-solution[0][i])<0.001);

  delete quad[0];
  quad[0] = NULL;

  MPI_Finalize();

  return 0;
}
