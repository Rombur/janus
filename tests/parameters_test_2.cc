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
  double inc_bottom(10.);
  BC_TYPE bottom_bc_type(most_normal);
  BC_TYPE right_bc_type(reflective);
  BC_TYPE top_bc_type(vacuum);
  BC_TYPE left_bc_type(reflective);
  FE_TYPE fe_type(pwld);
  QUAD_TYPE quad_type(ls);
  SOLVER_TYPE solver_type(gmres);
  d_vector src(3,0.);
  src[0] = 1.2;
  src[1] = 3.2;
  src[2] = 1.3;
  vector<d_vector> sigma_t(3,d_vector(3,0.));
  sigma_t[0][0] = 9.5;
  sigma_t[0][1] = 7;
  sigma_t[0][2] = 11.5;
  sigma_t[1][0] = 23.;
  sigma_t[1][1] = 18.;
  sigma_t[1][2] = 27.;
  sigma_t[2][0] = 34.5;
  sigma_t[2][1] = 27.;
  sigma_t[2][2] = 40.5;
  vector<vector<d_vector> > sigma_s(3,vector<d_vector>(3,d_vector(15,0.)));
  sigma_s[0][0][0] = 8.;
  sigma_s[0][0][1] = 7.;
  sigma_s[0][0][2] = 7.;
  sigma_s[0][0][3] = 5.;
  sigma_s[0][0][4] = 5.;
  sigma_s[0][0][5] = 5.;
  sigma_s[0][0][6] = 2.;
  sigma_s[0][0][7] = 2.;
  sigma_s[0][0][8] = 2.;
  sigma_s[0][0][9] = 2.;
  sigma_s[0][0][10] = -2.;
  sigma_s[0][0][11] = -2.;
  sigma_s[0][0][12] = -2.;
  sigma_s[0][0][13] = -2.;
  sigma_s[0][0][14] = -2.;
  sigma_s[0][1][0] = 5.5;
  sigma_s[0][1][1] = 4.5;
  sigma_s[0][1][2] = 4.5;
  sigma_s[0][1][3] = 2.5;
  sigma_s[0][1][4] = 2.5;
  sigma_s[0][1][5] = 2.5;
  sigma_s[0][2][0] = 10.;
  sigma_s[1][0][0] = 16.;
  sigma_s[1][0][1] = 14.;
  sigma_s[1][0][2] = 14.;
  sigma_s[1][0][3] = 10.;
  sigma_s[1][0][4] = 10.;
  sigma_s[1][0][5] = 10.;
  sigma_s[1][0][6] = 4.;
  sigma_s[1][0][7] = 4.;
  sigma_s[1][0][8] = 4.;
  sigma_s[1][0][9] = 4.;
  sigma_s[1][0][10] = -4.;
  sigma_s[1][0][11] = -4.;
  sigma_s[1][0][12] = -4.;
  sigma_s[1][0][13] = -4.;
  sigma_s[1][0][14] = -4.;
  sigma_s[1][1][0] = 11.;
  sigma_s[1][1][1] = 9.;
  sigma_s[1][1][2] = 9.;
  sigma_s[1][1][3] = 5.;
  sigma_s[1][1][4] = 5.;
  sigma_s[1][1][5] = 5.;
  sigma_s[1][2][0] = 20.;
  sigma_s[2][0][0] = 24.;
  sigma_s[2][0][1] = 21;
  sigma_s[2][0][2] = 21;
  sigma_s[2][0][3] = 15.;
  sigma_s[2][0][4] = 15.;
  sigma_s[2][0][5] = 15.;
  sigma_s[2][0][6] = 6.;
  sigma_s[2][0][7] = 6.;
  sigma_s[2][0][8] = 6.;
  sigma_s[2][0][9] = 6.;
  sigma_s[2][0][10] = -6.;
  sigma_s[2][0][11] = -6.;
  sigma_s[2][0][12] = -6.;
  sigma_s[2][0][13] = -6.;
  sigma_s[2][0][14] = -6.;
  sigma_s[2][1][0] = 16.5;
  sigma_s[2][1][1] = 13.5;
  sigma_s[2][1][2] = 13.5;
  sigma_s[2][1][3] = 7.5;
  sigma_s[2][1][4] = 7.5;
  sigma_s[2][1][5] = 7.5;
  sigma_s[2][2][0] = 30.;

  string filename("parameters_2.inp");

  PARAMETERS param(&filename);

  param.Read_parameters(3,3);

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
  assert(inc_bottom==param.Get_inc_bottom());

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
        j_max = 15;
      if (lvl==1)
        j_max = 6;
      if (lvl==2)
        j_max = 1;
      for (unsigned int j=0; j<j_max; ++j)
        assert(sigma_s[i][lvl][j]==param.Get_sigma_s(i)[lvl][j]);
    }

  return 0;
}
