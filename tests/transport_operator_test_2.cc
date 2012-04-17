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
#include "GLC.hh"
#include "LS.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
#include "TRANSPORT_OPERATOR.hh"
#include "TRIANGULATION.hh"

using namespace std;

typedef vector<double> d_vector;

int main(int argc,char** argv)
{
  unsigned int n_sources(2);
  unsigned int n_materials(2);
  unsigned int flux_moments_size(0);
  const double sqrt_4pi(2.*sqrt(M_PI));
  string geometry_inp("/home/bruno/Documents/Transport/janus/tests/geometry_transport_2.inp");
  string parameters_inp("/home/bruno/Documents/Transport/janus/tests/parameters_transport_2.inp");

  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  // Reference solution (scalar flux)
  d_vector solution(36,0.);
  solution[0] = -9.43310535e-2;
  solution[1] = 1.01687763e-1;
  solution[2] = 1.10930705e-1;
  solution[3] = 1.05328758e-1;
  solution[4] = -9.81478846e-2;
  solution[5] = 1.05558599e-1;
  solution[6] = 1.10930705e-1;
  solution[7] = 1.01687763e-1;
  solution[8] = -9.43310535e-2;
  solution[9] = 1.05558599e-1;
  solution[10] = -9.81478846e-2;
  solution[11] = 1.05328758e-1;
  solution[12] = 2.05948673;
  solution[13] = 2.05329349;
  solution[14] = 2.0587085;
  solution[15] = 2.0587085;
  solution[16] = 2.05329349;
  solution[17] = 2.05948673;
  solution[18] = 2.1154251;
  solution[19] = 2.1154251;

  TRIANGULATION triangulation(&geometry_inp);
  PARAMETERS parameters(&parameters_inp);
  vector<QUADRATURE*> quad;

  // Create the triangulation
  triangulation.Read_geometry();
  triangulation.Build_edges();

  // Create the parameters
  parameters.Read_parameters(n_sources,n_materials);

  // Build the quadrature
  unsigned int n_lvl(parameters.Get_n_levels());
  unsigned int tmp_sn(parameters.Get_sn_order());
  unsigned int tmp_L_max(parameters.Get_L_max());
  if (parameters.Get_mip()==true || parameters.Get_multigrid()==true)
    --n_lvl;
  quad.resize(n_lvl,NULL);
  for (unsigned int lvl=0; lvl<n_lvl; ++lvl)
  {
    if (parameters.Get_quad_type()==ls)
      quad[lvl] = new LS(tmp_sn,tmp_L_max,parameters.Get_galerkin());
    else
      quad[lvl] = new GLC(tmp_sn,tmp_L_max,parameters.Get_galerkin());
    quad[lvl]->Build_quadrature();

    tmp_sn = ceil(tmp_sn/2.);
    if (tmp_L_max>tmp_sn)
      tmp_L_max = tmp_sn;
  }                               

  // Build the dof handler
  DOF_HANDLER dof_handler(&triangulation,parameters);
  dof_handler.Compute_sweep_ordering(quad);
  
  // Flux moments map and vector
  flux_moments_size = dof_handler.Get_n_dof()*quad[0]->Get_n_mom();
  Epetra_Map flux_moments_map(flux_moments_size,0,comm);
  Epetra_MultiVector flux_moments(flux_moments_map,1);

  // Solve the transport
  const unsigned int lvl(0);
  const unsigned int max_lvl(parameters.Get_n_levels()-1);
  TRANSPORT_OPERATOR transport_operator(&dof_handler,&parameters,&quad,&comm,
      &flux_moments_map,lvl,max_lvl);

  // Compute right-hand side for GMRES
  Epetra_MultiVector rhs(flux_moments);
  transport_operator.Sweep(rhs,true);

  Epetra_LinearProblem problem(&transport_operator,&flux_moments,&rhs);

  AztecOO solver(problem);

  // Set Aztec parameters
  solver.SetAztecOption(AZ_precond,AZ_none);
  solver.SetAztecOption(AZ_solver,AZ_gmres);
  solver.SetAztecOption(AZ_conv,AZ_rhs);

  // Solve the transport equation
  solver.Iterate(parameters.Get_max_it(),parameters.Get_tolerance());

  // Apply the preconditioner to get the solution
  transport_operator.Apply_preconditioner(flux_moments);

  for (unsigned int i=0; i<20; ++i)
    assert(fabs(sqrt_4pi*flux_moments[0][i]-solution[i])<0.00001);


  MIP* mip(transport_operator.Get_mip());
  mip->Free_ml();
  
  for (unsigned int i=0; i<n_lvl; ++i)
  {
    delete quad[i];
    quad[i] = NULL;
  }

  MPI_Finalize();

  return 0;
}
