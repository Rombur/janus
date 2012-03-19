#include <cassert>
#include "PARAMETERS.hh"

using namespace std;

int main(int argc,char** argv)
{
  bool galerkin(false);
  bool multigrid(false);
  bool mip(false);
  unsigned int max_it(1000);
  unsigned int L_max(2);
  unsigned int sn(8);
  unsigned int n_levels(1);
  double tolerance(1e-6);
  double inc_bottom(-1.);
  double inc_right(0.);
  double inc_top(1.2);
  double inc_left(0.);
  FE_TYPE fe_type(bld);
  QUAD_TYPE quad_type(glc);
  SOLVER_TYPE solver_type(si);
  d_vector src(3,0.);
  src[0] = 0.;
  src[1] = 0.5;
  src[2] = 10.2;
  d_vector sigma_t(3,0.);
  sigma_t[0] = 10.1;
  sigma_t[1] = 11.5;
  sigma_t[2] = 12.5;
  vector<d_vector> sigma_s(3,d_vector(3,0.));
  sigma_s[0][0] = 10.;
  sigma_s[0][1] = 9.;
  sigma_s[0][2] = 8.;
  sigma_s[1][0] = 10.;
  sigma_s[1][1] = 9.; 
  sigma_s[1][2] = 8.; 
  sigma_s[2][0] = 10.;
  sigma_s[2][1] = 9.; 
  sigma_s[2][2] = 8.; 

  // Need to use AT_DATA instead of hard coding the path
  string filename("/home/bruno/Documents/Transport/janus/tests/parameters_1.txt");

  PARAMETERS param(&filename);

  param.Read_parameters(3,3);

  // Check mip flag
  assert(mip==param.Get_mip());

  // Check multigrid flag
  assert(multigrid==param.Get_multigrid());

  // Check the number of level
  assert(n_levels==param.Get_n_levels());

  // Check the incoming fluxes
  assert(inc_bottom==param.Get_inc_bottom());
  assert(inc_right==param.Get_inc_right());
  assert(inc_top==param.Get_inc_top());
  assert(inc_left==param.Get_inc_left());

  // Check the intensity of the source
  for (unsigned int i=0; i<3; ++i)
    assert(src[i]==param.Get_src(i));

  // Check the type of fe
  assert(fe_type==param.Get_fe_type());

  // Check the type of quadrature
  assert(quad_type==param.Get_quad_type());
  assert(galerkin==param.Get_galerkin());
  assert(L_max==param.Get_L_max());
  assert(sn==param.Get_sn_order());

  // Check the type of solver
  assert(solver_type==param.Get_solver_type());
  assert(max_it==param.Get_max_it());
  assert(tolerance==param.Get_tolerance());

  // Check the sigma_t
  for (unsigned int i=0; i<3; ++i)
    assert(sigma_t[i]==param.Get_sigma_t(i)[0]);

  // Check the sigma_s
  for (unsigned int i=0; i<3; ++i)
    for (unsigned int j=0; j<3; ++j)
      assert(sigma_s[i][j]==param.Get_sigma_s(i)[0][j]);

  return 0;
}
