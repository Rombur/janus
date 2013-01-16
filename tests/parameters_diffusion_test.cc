#include <cassert>
#include "PARAMETERS.hh"

using namespace std;

int main(int argc,char** argv)
{
  // mip_solver_type, inner_tolerance, group_tolerance, max_inner_it,
  // max_group_it, verbose, fe_type, src, incoming
  const unsigned int max_inner_it(100);
  const unsigned int max_group_it(10);
  const unsigned int n_groups(2);
  const unsigned int n_refinements(2);
  const unsigned int n_src(1);
  const double inner_tolerance(0.0001);
  const double group_tolerance(0.003);
  const double refinement_threshold(0.25);
  const BC_TYPE bottom_bc_type(vacuum);
  const BC_TYPE right_bc_type(isotropic);
  const BC_TYPE top_bc_type(vacuum);
  const BC_TYPE left_bc_type(isotropic);
  const FE_TYPE fe_type(bld);
  const MIP_SOLVER_TYPE mip_solver_type(agmg);
  d_vector inc_right(n_groups,0.);
  inc_right[0] = 2.;
  inc_right[1] = 3.;
  d_vector inc_left(n_groups,0.);
  inc_left[0] = 1.;
  inc_left[1] = 2.;
  vector<d_vector> src(n_src,d_vector (n_groups,0.));
  src[0][0] = 5.;
  src[0][1] = 10.;

  string filename("parameters_diffusion.inp");

  PARAMETERS param(&filename);

  param.Read_diffusion_parameters(n_src);

  // Check MIP solver
  assert(param.Get_mip_solver_type()==mip_solver_type);

  // Check inner tolerance
  assert(param.Get_inner_tolerance()==inner_tolerance);

  // Check group tolerance
  assert(param.Get_group_tolerance()==group_tolerance);

  // Check maximum number of inner iterations
  assert(param.Get_max_inner_it()==max_inner_it);

  // Check maximum number of group iterations
  assert(param.Get_max_group_it()==max_group_it);

  // Check the number of adaptive refinement to perform
  assert(param.Get_n_refinements()==n_refinements);

  // Check the refinement threshold
  assert(param.Get_refinement_threshold()==refinement_threshold);

  // Check the type of fe
  assert(fe_type==param.Get_fe_type());

  // Check the boundary condition types
  assert(bottom_bc_type==param.Get_bottom_bc_type());
  assert(right_bc_type==param.Get_right_bc_type());
  assert(top_bc_type==param.Get_top_bc_type());
  assert(left_bc_type==param.Get_left_bc_type());

  // Check the incoming fluxes
  for (unsigned int g=0; g<n_groups; ++g)
  {
    assert(inc_right[g]==param.Get_inc_right(g));
    assert(inc_left[g]==param.Get_inc_left(g));
  }

  return 0;
}
