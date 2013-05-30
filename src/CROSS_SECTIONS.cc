/*
Copyright (c) 2012, Bruno Turcksin.

This file is part of Janus.

Janus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
he Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Janus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Janus.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "CROSS_SECTIONS.hh"

CROSS_SECTIONS::CROSS_SECTIONS(string* cross_section_inputfile) :
  sigma_e_flag(false),
  L_max(0),
  n_g_groups(0),
  n_e_groups(0),
  n_p_groups(0),
  n_groups(0),
  n_supergroups(0),
  n_levels(0),
  cross_section_filename(cross_section_inputfile) 
{}

void CROSS_SECTIONS::Read_cepxs_cross_sections(const unsigned int n_mat, 
    const PERMUTATION_TYPE permutation_type)
{
  // Open the file to read it
  ifstream cross_section_file(cross_section_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(cross_section_file.good(),string("Unable to open the file ") +
        string(*cross_section_filename) + string(" containing the cross sections."));
  
  string skip_line;

  // Read the number of groups of gammas, electrons, positrons, and the number
  // of materials
  unsigned int n_materials(0);
  cross_section_file>>n_g_groups>>n_e_groups>>n_p_groups>>n_materials;
  n_groups = n_g_groups+n_e_groups+n_p_groups;
  
  // Skip the rest of line
  getline(cross_section_file,skip_line);

  // Check that n_mat is the same that n_materials
  Check(n_mat==n_materials,string("The number of materials in the geometry file ")
      +string(" and in the cross sections are different."));

  // Read L_max
  cross_section_file>>L_max;

  // Skip the rest of line
  getline(cross_section_file,skip_line);

  sigma_t.resize(n_mat,d_vector (n_groups,0.));
  sigma_e.resize(n_mat,d_vector (n_groups,0.));
  sigma_s.resize(n_mat,vector<vector<d_vector> > (n_groups, 
        vector<d_vector> (n_groups,d_vector (L_max+1,0.))));
  // Loop over the materials
  for (unsigned int i=0; i<n_mat; ++i)
  {
    // Store the sigma_t
    for (unsigned int g=0; g<n_groups; ++g)
      cross_section_file>>sigma_t[i][g];

    // Skip the rest of line
    getline(cross_section_file,skip_line);

    // Store the sigma_e
    for (unsigned int g=0; g<n_groups; ++g)
      cross_section_file>>sigma_e[i][g];

    // Skip the rest of line
    getline(cross_section_file,skip_line);

    // Store the sigma_s
    for (unsigned int mom=0; mom<L_max+1; ++mom)
      for (unsigned int g=0; g<n_groups; ++g)
        for (unsigned int gp=0; gp<n_groups; ++gp)
          cross_section_file>>sigma_s[i][g][gp][mom];
  }

  // Close the file
  cross_section_file.close();

  sigma_e_flag =true;

  // Apply the permutation on the cross sections
  Apply_cross_section_permutation(permutation_type,n_mat);
}

void CROSS_SECTIONS::Read_regular_cross_sections(const unsigned int n_mat,
    const PERMUTATION_TYPE permutation_type, const bool energy_deposition)
{
  // Open the file to read it
  ifstream cross_section_file(cross_section_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(cross_section_file.good(),string("Unable to open the file ") +
        string(*cross_section_filename) + string(" containing the cross sections."));
  
  // Read the number of groups, the number of materials, and L_max
  unsigned int n_materials(0);
  cross_section_file>>n_groups>>n_materials>>L_max;

  // Check that n_mat is the same that n_materials
  Check(n_mat==n_materials,string("The number of materials in the geometry file ")
      +string(" and in the cross sections are different."));

  sigma_t.resize(n_mat,d_vector (n_groups,0.));
  sigma_s.resize(n_mat,vector<vector<d_vector> > (n_groups, 
        vector<d_vector> (n_groups,d_vector (L_max+1,0.))));
  // Loop over the materials
  for (unsigned int i=0; i<n_mat; ++i)
  {
    // Store the sigma_t
    for (unsigned int g=0; g<n_groups; ++g)
      cross_section_file>>sigma_t[i][g];

    // Store the sigma_e if present
    if (energy_deposition==true)
    {
      sigma_e.resize(n_mat,d_vector (n_groups,0.));
      for (unsigned int g=0; g<n_groups; ++g)
        cross_section_file>>sigma_e[i][g];
    }

    // Store the sigma_s
    for (unsigned int g=0; g<n_groups; ++g)
      for (unsigned int gp=0; gp<n_groups; ++gp)
        for (unsigned int mom=0; mom<L_max+1; ++mom)
          cross_section_file>>sigma_s[i][g][gp][mom];
  }

  // Close the file
  cross_section_file.close();

  sigma_e_flag = energy_deposition;

  // Apply the permutation on the cross sections
  Apply_cross_section_permutation(permutation_type,n_mat);
}

void CROSS_SECTIONS::Build_fokker_planck_xs(const unsigned int n_mat)
{
  // Open the file to read it
  ifstream cross_section_file(cross_section_filename->c_str(),ios::in);
  

  // Check that the file was open properly
  Check(cross_section_file.good(),string("Unable to open the file ") +
        string(*cross_section_filename) + string(" containing the cross sections."));
  
  // Read the number of materials and L_max
  unsigned int n_materials(0);
  cross_section_file>>n_materials>>L_max;

  // When Fokker-Planck cross sections are used, only one group is allowed
  n_groups = 1;
  n_supergroups = 1;

  // Check that n_mat is the same that n_materials
  Check(n_mat==n_materials,string("The number of materials in the geometry file ")
      +string(" and in the cross sections are different."));

  d_vector alpha(n_mat,0.);
  sigma_t.resize(n_mat,d_vector (1,0.));
  sigma_s.resize(n_mat,vector<vector<d_vector> > (1,vector<d_vector> (1,
          d_vector (L_max+1,0.))));
  // Loop over the materials
  for (unsigned int i=0; i<n_mat; ++i)
  {
    // Store the sigma_t
    cross_section_file>>sigma_t[i][0];

    // Store the alpha
    cross_section_file>>alpha[i];
  }

  // Close the file
  cross_section_file.close();

  // Build the Fokker-Planck cross sections 
  for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
    for (unsigned int l=0; l<L_max+1; ++l)
      sigma_s[i_mat][0][0][l] = alpha[i_mat]/2.*(L_max*(L_max+1)-l*(l+1));
}

void CROSS_SECTIONS::Apply_ang_lvls_and_tc(const bool ang_lvls,const bool tc,
    const bool optimal,const unsigned int n_mat,const unsigned int sn)
{  
  if (ang_lvls==true)
  {
    Check(ang_lvls==tc,"Angular multigrid has to be used with transport correction.");
    Apply_angular_levels(n_mat,sn);
  }
  else
  {
    n_levels = 1;
    sigma_t_lvl.resize(n_mat,vector<d_vector> (n_groups,d_vector(n_levels,0.)));
    sigma_s_lvl.resize(n_mat,vector<vector<vector<d_vector> > >(n_groups,
           vector<vector<d_vector> > (n_groups,vector<d_vector>(n_levels,
               d_vector(L_max+1,0.)))));
    for (unsigned int material_id=0; material_id<n_mat; ++material_id)
      for (unsigned int g=0; g<n_groups; ++g)
        sigma_t_lvl[material_id][g][0] = sigma_t[material_id][g];
    for (unsigned int material_id=0; material_id<n_mat; ++material_id)
      for (unsigned int g=0; g<n_groups; ++g)
        for (unsigned int gp=0; gp<n_groups; ++gp)
          for (unsigned int l=0; l<=L_max; ++l)
            sigma_s_lvl[material_id][g][gp][0][l] = sigma_s[material_id][g][gp][l];
  }
  
  if (tc==true)
  {
    // "standard" extended correction is used, the last sigma_s,l is
    // used as correction and L_max is decreased by one.
    if (optimal==false)
      --L_max;
    
    for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
    {
      for (unsigned int g=0; g<n_groups; ++g)
      {
        double correction(0.);
        double L(L_max); 
        for (unsigned int lvl=0; lvl<n_levels; ++lvl)
        {
          double L2(ceil(L/2.));
          if (optimal==true)
          {
            if (L==0.)
              correction = 0.;
            else
            {
              correction = (sigma_s_lvl[i_mat][g][g][lvl][L]+
                  sigma_s_lvl[i_mat][g][g][lvl][L2+1])/2.;
            }
          }
          else
          {
            if (lvl==0)
              correction = sigma_s[i_mat][g][g][L_max+1];
            else
              correction = sigma_s[i_mat][g][g][L+1]; 
          }
          Apply_transport_correction(i_mat,g,lvl,correction);
          L = L2;
        }
      }
    }
  }
}

void CROSS_SECTIONS::Apply_angular_levels(const unsigned int n_mat,
    const unsigned int sn)
{
  n_levels = ceil(log(double(sn))/log(2.));
  // Loop over the level of the angular multigrid
  sigma_t_lvl.resize(n_mat,vector<d_vector> (n_groups,d_vector(n_levels,0.)));
  sigma_s_lvl.resize(n_mat,vector<vector<vector<d_vector> > >(n_groups,
        vector<vector<d_vector> > (n_groups,vector<d_vector>(n_levels,
            d_vector(L_max+1,0.)))));
  for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
  {
    for (unsigned int g=0; g<n_groups; ++g)
    {
      for (unsigned int lvl=0; lvl<n_levels; ++lvl)
      {
        sigma_t_lvl[i_mat][g][lvl] = sigma_t[i_mat][g];
        for (unsigned int gp=0; gp<n_groups; ++gp)
          sigma_s_lvl[i_mat][g][gp][lvl] = d_vector(sigma_s[i_mat][g][gp].begin(),
              sigma_s[i_mat][g][gp].end());
      }
    }
  }
}

void CROSS_SECTIONS::Apply_transport_correction(const unsigned int i_mat,
    const unsigned int g,const unsigned int lvl,const double correction)
{
  sigma_t_lvl[i_mat][g][lvl] -= correction;
  d_vector::iterator sigma_s(sigma_s_lvl[i_mat][g][g][lvl].begin());
  d_vector::iterator sigma_s_end(sigma_s_lvl[i_mat][g][g][lvl].end());
  for (; sigma_s<sigma_s_end; ++sigma_s)
    *sigma_s -= correction;
}

void CROSS_SECTIONS::Apply_cross_section_permutation(
    const PERMUTATION_TYPE permutation_type,const unsigned int n_mat)
{
  bool identity(true);
  ui_vector permutation_v(n_groups,0.);
  vector<d_vector> permutation_m(n_groups,d_vector (n_groups,0.));

  switch (permutation_type)
  {
    case (linear) :
      {
        Create_linear_permutation(permutation_m,permutation_v,identity);
        break;
      }
    case (logarithmic) :
      {
        Create_log_permutation(permutation_m,permutation_v);
        identity = false;
        break;
      }
    case (none) :
      {
        n_supergroups = n_groups;
        break;
      }
    default :
      {
        Check(false,"Unknown permutation type.");
      }
  }
  if (!identity)
  {
    Matrix_multiplication(n_mat,permutation_m);
    Vector_permutation(n_mat,permutation_v);
  }
}

unsigned int CROSS_SECTIONS::Compute_gcd(unsigned int n_a,unsigned int n_b)
{   
  unsigned int tmp;
  while (n_b)
  {
    tmp = n_b;
    n_b = n_a%n_b;
    n_a = tmp;
  }

  return n_a;
}

void CROSS_SECTIONS::Matrix_multiplication(const unsigned int n_mat,
    vector<d_vector> const &permutation_m)
{
  vector<vector<vector<d_vector> > > result(n_mat,vector<vector<d_vector> >
      (n_groups,vector<d_vector> (n_groups,d_vector(L_max+1,0.))));

  // result = XS * P
  for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
    for (unsigned int g=0; g<n_groups; ++g)
      for (unsigned int gp=0; gp<n_groups; ++gp)
        for (unsigned int gpp=0; gpp<n_groups; ++gpp)
          for (unsigned int l=0; l<L_max+1; ++l)
            result[i_mat][g][gpp][l] += sigma_s[i_mat][g][gp][l]*
              permutation_m[gp][gpp];

  sigma_s.clear();
  sigma_s.resize(n_mat,vector<vector<d_vector> > (n_groups,
        vector<d_vector> (n_groups,d_vector(L_max+1,0.))));

  // XS = P^t * result
  for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
    for (unsigned int g=0; g<n_groups; ++g)
      for (unsigned int gp=0; gp<n_groups; ++gp)
        for (unsigned int gpp=0; gpp<n_groups; ++gpp)
          for (unsigned int l=0; l<L_max+1; ++l)
            sigma_s[i_mat][g][gpp][l] += permutation_m[gp][g]*
              result[i_mat][gp][gpp][l];
}

void CROSS_SECTIONS::Vector_permutation(const unsigned int n_mat,
    ui_vector const &permutation_v)
{
  for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
  {
    d_vector e_tmp(sigma_e[i_mat]);
    d_vector t_tmp(sigma_t[i_mat]);
    for (unsigned int g=0; g<n_g_groups; ++g)
    {
      sigma_e[i_mat][permutation_v[g]] = e_tmp[g];
      sigma_t[i_mat][permutation_v[g]] = t_tmp[g];
    }
  }
}

void CROSS_SECTIONS::Create_linear_permutation(vector<d_vector> &permutation_m,
    ui_vector &permutation_v,bool &identity)
{
  const unsigned int gcd(Compute_gcd(n_g_groups,n_e_groups));
  n_supergroups = n_groups/gcd;
  if (gcd!=1)
  {
    identity = false;
    const unsigned int gcd_g(n_g_groups/gcd);
    const unsigned int gcd_e(n_e_groups/gcd);
    
    for (unsigned int i=0; i<n_g_groups; ++i)
    {
      unsigned int column(i/gcd_g*(gcd_g+gcd_e)+i%gcd_g);
      permutation_m[i][column] = 1.0;
      permutation_v[i] = column;
    }

    for (unsigned int i=0; i<n_e_groups; ++i)
    {
      unsigned int column(gcd_g+i/gcd_e*(gcd_g+gcd_e)+i%gcd_e);
      permutation_m[i+n_g_groups][column] = 1.0;
      permutation_v[i+n_g_groups] = column;
    }
  }
}

void CROSS_SECTIONS::Create_log_permutation(vector<d_vector> &permutation_m,
    ui_vector &permutation_v)
{
  for (unsigned int i=0; i<n_g_groups; i++)
  {
    if (i<n_e_groups)
    {
      permutation_m[i][2*i] = 1.0;
      permutation_v[i] = 2*i;
    }
    else
    {
      permutation_m[i][2*n_e_groups+i] = 1.0;
      permutation_v[i] = 2*n_e_groups+i;
    }
  }

  for (unsigned int i=0; i<n_e_groups; i++)
  {
    if (i<n_g_groups)
    {
      permutation_m[i+n_g_groups][2*i+1] = 1.0;
      permutation_v[i+n_g_groups] = 2*i+1;
    }
    else
    {
      permutation_m[i+n_g_groups][n_g_groups+i] = 1.0;
      permutation_v[i+n_g_groups] = n_g_groups+i;
    }
  }
}
