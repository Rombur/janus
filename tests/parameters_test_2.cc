#include <cassert>
#include "PARAMETERS.hh"

using namespace std;

int main(int argc,char** argv)
{
  const bool galerkin(true);
  const bool multigrid(true);
  const bool mip(true);
  const unsigned int max_inner_it(270);
  const unsigned int max_supergroup_it(20);
  const unsigned int max_group_it(2);
  const unsigned int sn(4);
  const unsigned int n_levels(2);
  const unsigned int n_groups(3);
  const unsigned int n_src(3);
  const double inner_tolerance(1e-10);
  const double group_tolerance(1e-8);
  const double weight_sum(2.*M_PI);
  const BC_TYPE bottom_bc_type(most_normal);
  const BC_TYPE right_bc_type(reflective);
  const BC_TYPE top_bc_type(vacuum);
  const BC_TYPE left_bc_type(reflective);
  const FE_TYPE fe_type(pwld);
  const QUAD_TYPE quad_type(ls);
  const SOLVER_TYPE solver_type(gmres);
  d_vector inc_bottom(n_groups,0.);
  inc_bottom[0] = 10.;
  inc_bottom[1] = 5.;
  inc_bottom[2] = 2.5;
  vector<d_vector> src(n_src,d_vector (n_groups,0.));
  src[0][0] = 1.2;
  src[0][1] = 3.2;
  src[0][2] = 1.3;
  src[1][0] = 3.2;
  src[1][1] = 4.2;
  src[1][2] = 5.2;
  src[2][0] = 1.3;
  src[2][1] = 3.3;
  src[2][2] = 2.3;

  string filename("parameters_2.inp");

  PARAMETERS param(&filename);

  param.Read_parameters(n_src);

  // Check sum of the weights of the quadrature
  assert((weight_sum-param.Get_weight_sum())<1e-14);

  // Check mip flag
  assert(mip==param.Get_mip());

  // Check multigrid flag
  assert(multigrid==param.Get_multigrid());

  // Check the number of level
  assert(n_levels==param.Get_n_levels());

  // Check the boundary condition types
  assert(bottom_bc_type==param.Get_bottom_bc_type());
  assert(right_bc_type==param.Get_right_bc_type());
  assert(top_bc_type==param.Get_top_bc_type());
  assert(left_bc_type==param.Get_left_bc_type());
  
  // Check the incoming flux
  for (unsigned int g=0; g<n_groups; ++g)
    assert(inc_bottom[g]==param.Get_inc_bottom(g));

  // Check the intensity of the source
  for (unsigned int i=0; i<n_src; ++i)
  {
    d_vector source(param.Get_src(i));
    for (unsigned int g=0; g<n_groups; ++g)
      assert(src[i][g]==source[g]);
  }

  // Check the type of fe
  assert(fe_type==param.Get_fe_type());

  // Check the type of quadrature
  assert(quad_type==param.Get_quad_type());
  assert(galerkin==param.Get_galerkin());
  assert(sn==param.Get_sn_order());

  // Check the type of solver
  assert(solver_type==param.Get_solver_type());
  assert(max_inner_it==param.Get_max_inner_it());
  assert(max_supergroup_it==param.Get_max_supergroup_it());
  assert(max_group_it==param.Get_max_group_it());
  assert(inner_tolerance==param.Get_inner_tolerance());
  assert(group_tolerance==param.Get_group_tolerance());

  return 0;
}
