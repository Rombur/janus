#include <cassert>
#include <string>
#include "CROSS_SECTIONS.hh"

using namespace std;

typedef vector<double> d_vector;

int main()
{
  const bool ang_lvls(false);
  const bool tc(false);
  const bool optimal(false);
  const unsigned int n_levels(1);
  const unsigned int sn(8);
  const unsigned int n_mat(2);
  const unsigned int L_max(8);
  const unsigned int n_groups(1);
  const unsigned int n_supergroups(1);
  const unsigned int n_grps_in_supergrp(1);
  string cross_sections_inp("cross_sections_fp.inp");
  vector<vector<d_vector> > sigma_t_ref(n_mat,vector<d_vector> (n_groups,
        d_vector(n_levels,0.)));
  sigma_t_ref[0][0][0] = 36.;
  sigma_t_ref[1][0][0] = 20.;
  vector<vector<vector<vector<d_vector> > > > sigma_s_ref(n_mat,
      vector<vector<vector<d_vector> > > (n_groups,
        vector<vector<d_vector> > (n_groups,vector<d_vector> (n_levels,
            d_vector(L_max+1,0.)))));
  sigma_s_ref[0][0][0][0][0] = 36.;
  sigma_s_ref[0][0][0][0][1] = 35.;
  sigma_s_ref[0][0][0][0][2] = 33.;
  sigma_s_ref[0][0][0][0][3] = 30.;
  sigma_s_ref[0][0][0][0][4] = 26.;
  sigma_s_ref[0][0][0][0][5] = 21.;
  sigma_s_ref[0][0][0][0][6] = 15.;
  sigma_s_ref[0][0][0][0][7] = 8.;
  sigma_s_ref[1][0][0][0][0] = 18.;
  sigma_s_ref[1][0][0][0][1] = 17.5;
  sigma_s_ref[1][0][0][0][2] = 16.5;
  sigma_s_ref[1][0][0][0][3] = 15.;
  sigma_s_ref[1][0][0][0][4] = 13.;
  sigma_s_ref[1][0][0][0][5] = 10.5;
  sigma_s_ref[1][0][0][0][6] = 7.5;
  sigma_s_ref[1][0][0][0][7] = 4.;
  CROSS_SECTIONS cross_sections(&cross_sections_inp);

  // Build the Fokker-Planck cross sections.
  cross_sections.Build_fokker_planck_xs(n_mat);

  // Check L max
  assert(L_max==cross_sections.Get_L_max());

  // Check the number of groups
  assert(n_groups==cross_sections.Get_n_groups());

  // Check the number of supergroups
  assert(n_supergroups==cross_sections.Get_n_supergroups());

  // Check the number of groups in a supergroup
  assert(n_grps_in_supergrp==cross_sections.Get_n_grps_in_supergrp());

  // Build sigma_s_lvl and sigma_t_lvl when no multigrid or extended transport
  // correction is used
  cross_sections.Apply_ang_lvls_and_tc(ang_lvls,tc,optimal,n_mat,sn);

  // Check the number of levels
  assert(n_levels==cross_sections.Get_n_levels());
  
  // Check the total cross sections and the scattering cross sections
  for (unsigned int material_id=0; material_id<n_mat; ++material_id)
  {
    vector<d_vector> sigma_t(cross_sections.Get_sigma_t(material_id));
    vector<vector<vector<d_vector> > > sigma_s(cross_sections.Get_sigma_s(
          material_id));
    for (unsigned int g=0; g<n_groups; ++g)
      for (unsigned int lvl=0; lvl<n_levels; ++lvl)
        assert(sigma_t[g][lvl]==sigma_t_ref[material_id][g][lvl]);
    for (unsigned int g=0; g<n_groups; ++g)
      for (unsigned int gp=0; gp<n_groups; ++gp)
        for (unsigned int lvl=0; lvl<n_levels; ++lvl)
          for (unsigned int l=0; l<=L_max; ++l)
            assert(sigma_s[g][gp][lvl][l]==sigma_s_ref[material_id][g][gp][lvl][l]);
  }
  
  return 0;
}
