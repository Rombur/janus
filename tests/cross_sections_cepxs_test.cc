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
  const unsigned int sn(2);
  const unsigned int n_mat(1);
  const unsigned int L_max(2);
  const unsigned int n_groups(2);
  const unsigned int n_supergroups(2);
  const unsigned int n_grps_in_supergrp(1);
  string cross_sections_inp("cross_sections_cepxs.inp");
  vector<vector<d_vector> > sigma_t_ref(n_mat,vector<d_vector> (n_groups,
        d_vector(n_levels,0.)));
  sigma_t_ref[0][0][0] = 2.;
  sigma_t_ref[0][1][0] = 3.;
  vector<d_vector> sigma_e_ref(n_mat,d_vector (n_groups,0.));
  sigma_e_ref[0][0] = 5.;
  sigma_e_ref[0][1] = 7.;
  vector<vector<vector<vector<d_vector> > > > sigma_s_ref(n_mat,
      vector<vector<vector<d_vector> > > (n_groups,
        vector<vector<d_vector> > (n_groups,vector<d_vector> (n_levels,
            d_vector(L_max+1,0.)))));
  sigma_s_ref[0][0][0][0][0] = 23.;
  sigma_s_ref[0][0][0][0][1] = 29.;
  sigma_s_ref[0][0][0][0][2] = 31.;
  sigma_s_ref[0][0][1][0][0] = 37.;
  sigma_s_ref[0][0][1][0][1] = 41.;
  sigma_s_ref[0][0][1][0][2] = 43.;
  sigma_s_ref[0][1][0][0][0] = 47.;
  sigma_s_ref[0][1][0][0][1] = 53.;
  sigma_s_ref[0][1][0][0][2] = 59.;
  sigma_s_ref[0][1][1][0][0] = 61.;
  sigma_s_ref[0][1][1][0][1] = 67.;
  sigma_s_ref[0][1][1][0][2] = 71.;

  CROSS_SECTIONS cross_sections(&cross_sections_inp);

  // Read the cross sections.
  cross_sections.Read_cepxs_cross_sections(n_mat,none);

  // Check the number of groups
  assert(n_groups==cross_sections.Get_n_groups());

  // Check the number of supergroups
  assert(n_supergroups==cross_sections.Get_n_supergroups());

  // Check the number of groups in a supergroup
  assert(n_grps_in_supergrp==cross_sections.Get_n_grps_in_supergrp());

  // Build sigma_s_lvl and sigma_t_lvl when no multigrid or extended transport
  // correction is used
  cross_sections.Apply_ang_lvls_and_tc(ang_lvls,tc,optimal,n_mat,sn);

  // Check L max
  assert(L_max==cross_sections.Get_L_max());

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
