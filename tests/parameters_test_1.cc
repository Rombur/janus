#include <cassert>
#include "PARAMETERS.hh"

using namespace std;

typedef vector<double> d_vector;

int main(int argc,char** argv)
{
  const bool galerkin(false);
  const bool multigrid(false);
  const bool mip(false);
  const unsigned int max_inner_it(1000);
  const unsigned int max_supergroup_it(10);
  const unsigned int max_group_it(20);
  const unsigned int sn(8);
  const unsigned int n_levels(1);
  const unsigned int n_groups(2);
  const unsigned int n_src(2);
  const unsigned int verbose(0);
  const double inner_tolerance(1e-8);
  const double group_tolerance(1e-6);
  const double weight_sum(4.*M_PI);
  const BC_TYPE bottom_bc_type(most_normal);
  const BC_TYPE right_bc_type(isotropic);
  const BC_TYPE top_bc_type(isotropic);
  const BC_TYPE left_bc_type(most_normal);
  const FE_TYPE fe_type(bld);
  const QUAD_TYPE quad_type(glc);
  const SOLVER_TYPE solver_type(si);
  const XS_TYPE xs_type(fp);
  d_vector inc_bottom(n_groups,0.);
  inc_bottom[0] = 1.0;
  inc_bottom[1] = 0.5;
  inc_bottom[2] = 0.25;
  d_vector inc_right(n_groups,0.);
  inc_right[0] = 0.3;
  inc_right[1] = 0.15;
  inc_right[2] = 0.075;
  d_vector inc_top(n_groups,0.);
  inc_top[0] = 1.2;
  inc_top[1] = 0.6;
  inc_top[2] = 0.3;
  d_vector inc_left(n_groups,0.);
  inc_left[0] = 0.4;
  inc_left[1] = 0.2;
  inc_left[2] = 0.1;
  vector<d_vector> src(n_src,d_vector (n_groups,0.));
  src[0][0] = 0.;
  src[0][1] = 0.5;
  src[0][2] = 10.2;
  src[1][0] = 2.3;
  src[1][1] = 4.0;
  src[1][2] = 6.9;

  string filename("parameters_1.inp");

  PARAMETERS param(&filename);

  param.Read_parameters(n_src);

  // Check sum of the weights of the quadrature
  assert((weight_sum-param.Get_weight_sum())<1e-14);

  // Check the type of cross section
  assert(xs_type==param.Get_xs_type());

  // Check verbose flag
  assert(verbose==param.Get_verbose());

  // Check mip flag
  assert(mip==param.Get_mip());

  // Check multigrid flag
  assert(multigrid==param.Get_multigrid());

  // Check the type of quadrature
  assert(quad_type==param.Get_quad_type());
  assert(galerkin==param.Get_galerkin());
  assert(sn==param.Get_sn_order());

  // Check the type of fe
  assert(fe_type==param.Get_fe_type());

  // Check the type of solver
  assert(solver_type==param.Get_solver_type());
  assert(max_inner_it==param.Get_max_inner_it());
  assert(max_supergroup_it==param.Get_max_supergroup_it());
  assert(max_group_it==param.Get_max_group_it());
  assert(inner_tolerance==param.Get_inner_tolerance());

  // Check the number of level
  assert(n_levels==param.Get_n_levels());

  // Check the intensity of the source
  for (unsigned int i=0; i<n_src; ++i)
  {
    d_vector source(param.Get_src(i));
    for (unsigned int g=0; g<n_groups; ++g)
      assert(src[i][g]==source[g]);
  }
  assert(group_tolerance==param.Get_group_tolerance());

  // Check the boundary condition types
  assert(bottom_bc_type==param.Get_bottom_bc_type());
  assert(right_bc_type==param.Get_right_bc_type());
  assert(top_bc_type==param.Get_top_bc_type());
  assert(left_bc_type==param.Get_left_bc_type());

  // Check the incoming fluxes
  for (unsigned int g=0; g<n_groups; ++g)
  {
    assert(inc_bottom[g]==param.Get_inc_bottom(g));
    assert(inc_right[g]==param.Get_inc_right(g));
    assert(inc_top[g]==param.Get_inc_top(g));
    assert(inc_left[g]==param.Get_inc_left(g));
  }

  return 0;
}
