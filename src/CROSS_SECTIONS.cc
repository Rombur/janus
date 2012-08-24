#include "CROSS_SECTIONS.cc"

CROSS_SECTIONS::CROSS_SECTIONS(string* cross_section_inputfile) :
  n_g_groups(0),
  n_e_groups(0),
  n_p_groups(0),
  cross_section_filename(cross_section_inputfile) 
{}

void CROSS_SECTIONS::Read_cepxs_cross_sections(const unsigned int n_mat)
{
  // Open the file to read it
  ifstream cross_section_file(cross_section_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(cross_section_file,good(),string("Unable to open the file " +
        *cross_section_filename + " containing the cross sections."));
  
  string skip_line;

  // Read the number of groups of gammas, electrons, positrons, and the number
  // of materials
  unsigned int n_materials(0);
  cross_section_file>>n_g_groups>>n_e_groups>>n_p_groups>>n_materials;
  n_groups = n_g_groups+n_e_groups+n_p_groups;
  
  // Skip the rest of line
  getline(cross_section_file,skip_line);

  // Check that n_mat is the same that n_materials
  Check(n_mat==n_materials,string("The number of materials in the geometry file "+
        " and in the cross sections are different."));

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
          cross_section_file>>sigma_s[g][gp][mom];
  }

  // Close the file
  cross_section_file.close();
}

void CROSS_SECTIONS::Read_regular_cross_sections(const unsigned int n_mat,
    bool energy_deposition)
{
  // Open the file to read it
  ifstream cross_section_file(cross_section_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(cross_section_file,good(),string("Unable to open the file " +
        *cross_section_filename + " containing the cross sections."));
  
  // Read the number of groups, the number of materials, and L_max
  unsigned int n_materials(0);
  cross_section_file>>n_groups>>n_materials>>L_max;

  // Check that n_mat is the same that n_materials
  Check(n_mat==n_materials,string("The number of materials in the geometry file "+
        " and in the cross sections are different."));

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
    for (unsigned int mom=0; mom<L_max+1; ++mom)
      for (unsigned int g=0; g<n_groups; ++g)
        for (unsigned int gp=0; gp<n_groups; ++gp)
          cross_section_file>>sigma_s[g][gp][mom];
  }

  // Close the file
  cross_section_file.close();
}

void CROSS_SECTIONS::Build_fokker_planck_xs(const unsigned int n_mat)
{
  // Open the file to read it
  ifstream cross_section_file(cross_section_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(cross_section_file,good(),string("Unable to open the file " +
        *cross_section_filename + " containing the cross sections."));
  
  // Read the number of materials
  unsigned int n_materials(0);
  cross_section_file>>n_materials;

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
  {
    unsigned int index(0);
    for (unsigned int l=0; l<L_max+1; ++l)
        sigma_s[i_mat][index] = alpha[i_mat]/2.*(L_max*(L_max+1)-l*(l+1));
  }
}

void CROSS_SECTIONS::Apply_ang_lvls_and_tc(const bool ang_lvls,const bool tc,
    const bool optimal)
{  
  if (ang_lvls==true)
  {
    Check(ang_lvls==tc,"Angular multigrid has to be used with "+
        "transport correction.");
    Apply_angular_levels();
  }
  else
    n_levels = 1;
  
  if (tc==true)
  {
    // e "standard" extended correction is used, the last sigma_s,l is
    // used as correction and L_max is decreased by one.
    if (optimal==false)
      --L_max;
    
    for (unsigned int i_mat=0; i_mat<n_mat; ++i)
    {
      for (unsigned int g=0; g<n_groups; ++g)
      {
        double L(L_max); 
        for (unsigned int lvl=0; lvl<n_levels; ++lvl)
        {
          double L_2(ceil(L/2.));
          if (optimal==false)
          {
            if (lvl==0)
              correction = sigma_s[i_mat][g][g][L_max+1];
            else
              correction = sigma_s_lvl[i_mat][g][g][lvl-1][L_2]; 
          }
        }
        Apply_transport_correction(i_mat,g,lvl,L2,correction)
      }
      L = ceil(L/2.);
      if (L==1.)
        L = 0.;
    }
  }
}

void CROSS_SECTIONS::Apply_angular_levels(sn)
{
  n_levels = ceil(log(double(sn))/log(2.));
  // Loop over the level of the angular multigrid
  sigma_t_lvl.resize(n_mat,vector<d_vector> (n_groups,d_vector(n_levels,0.)));
  sigma_s_lvl.resize(n_mat,vector<vector<d_vector> >(n_groups,
        vector<d_vector> (n_groups,d_vector(n_levels,0.))));
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

void CROSS_SECTIONS::Apply_transport_correction(unsigned int i_mat,unsigned int g,
    unsigned int lvl, double L,double correction)
{
  if (optimal==true)
  {
    if (L==0.)
      correction = 0.;
    else
    {
      if (L==1.)
        L = 0.;
      unsigned int pos(0);
      for (unsigned int jj=0; jj<=L; ++jj)
        pos += jj+1;
      correction = (sigma_s_lvl[i_mat][g][g][lvl][pos]+
          sigma_s_lvl[i_mat][g][g][lvl].back())/2.;
    }
  }
  
  sigma_t_lvl[i_mat][g][lvl] -= correction;
  for (unsigned int gp=0; gp<n_groups; ++gp)
  {
    d_vector::iterator sigma_s(sigma_s_lvl[i_mat][g][gp][lvl].begin());
    d_vector::iterator sigma_s_end(sigma_s_lvl[i_mat][g][gp][lvl].end());
    for (; sigma_s<sigma_s_end; ++sigma_s)
      *sigma_s -= correction;
  }
}
