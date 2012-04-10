#include <cassert>
#include "PARAMETERS.hh"

using namespace std;

int main(int argc,char** argv)
{
  bool galerkin(true);
  bool multigrid(true);
  bool mip(true);
  unsigned int max_it(270);
  unsigned int L_max(4);
  unsigned int sn(4);
  unsigned int n_levels(3);
  double tolerance(1e-10);
  double inc_bottom(-1.);
  double inc_right(0.);
  double inc_top(-1.);
  double inc_left(0.);
  FE_TYPE fe_type(pwld);
  QUAD_TYPE quad_type(ls);
  SOLVER_TYPE solver_type(gmres);
  d_vector src(3,0.);
  src[0] = 1.2;
  src[1] = 3.2;
  src[2] = 1.3;
  vector<d_vector> sigma_t(3,d_vector(3,0.));
  sigma_t[0][0] = 8.;
  sigma_t[0][1] = 3.5;
  sigma_t[0][2] = 11.5;
  sigma_t[1][0] = 20.;
  sigma_t[1][1] = 11.;
  sigma_t[1][2] = 27.;
  sigma_t[2][0] = 30.;
  sigma_t[2][1] = 16.5;
  sigma_t[2][2] = 40.5;
  vector<vector<d_vector> > sigma_s(3,vector<d_vector>(3,d_vector(5,0.)));
  sigma_s[0][0][0] = 6.5;
  sigma_s[0][0][1] = 5.5;
  sigma_s[0][0][2] = 3.5;
  sigma_s[0][0][3] = 0.5;
  sigma_s[0][0][4] = -3.5;
  sigma_s[0][1][0] = 2.;
  sigma_s[0][1][1] = 1.;
  sigma_s[0][1][2] = -1.;
  sigma_s[0][2][0] = 10.;
  sigma_s[1][0][0] = 13.;
  sigma_s[1][0][1] = 11.;
  sigma_s[1][0][2] = 7.;
  sigma_s[1][0][3] = 1.;
  sigma_s[1][0][4] = -7.;
  sigma_s[1][1][0] = 4.;
  sigma_s[1][1][1] = 2.;
  sigma_s[1][1][2] = -2.;
  sigma_s[1][2][0] = 20.;
  sigma_s[2][0][0] = 19.5;
  sigma_s[2][0][1] = 16.5;
  sigma_s[2][0][2] = 10.5;
  sigma_s[2][0][3] = 1.5;
  sigma_s[2][0][4] = -10.5;
  sigma_s[2][1][0] = 6.;
  sigma_s[2][1][1] = 3.;
  sigma_s[2][1][2] = -3.;
  sigma_s[2][2][0] = 30.;

  // Need to use AT_DATA instead of hard coding the path
  string filename("/home/bruno/Documents/Transport/janus/tests/parameters_2.inp");

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
    for (unsigned int lvl=0; lvl<3; ++lvl)
      assert(sigma_t[i][lvl]==param.Get_sigma_t(i)[lvl]);

  // Check the sigma_s
  for (unsigned int i=0; i<3; ++i)
    for (unsigned int lvl=0; lvl<3; ++lvl)
    {
      unsigned int j_max(0);
      if (lvl==0)
        j_max = 5;
      if (lvl==1)
        j_max = 2;
      if (lvl==2)
        j_max = 1;
      for (unsigned int j=0; j<j_max; ++j)
        assert(sigma_s[i][lvl][j]==param.Get_sigma_s(i)[lvl][j]);
    }

  return 0;
}
