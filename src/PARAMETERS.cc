#include "PARAMETERS.hh"

PARAMETERS::PARAMETERS(string* parameters_inputfile) :
  parameters_filename(parameters_inputfile)
{}

void PARAMETERS::Read_parameters(unsigned int n_src,unsigned int n_mat)
{
  // Open the file to read it
  ifstream parameters_file(parameters_filename->c_str(),ios::in);

  string fe_type_str;
  string fokker_planck_str;
  string galerkin_str;
  string mip_str;
  string mip_solver_type_str;
  string multigrid_str;
  string optimal_str;
  string quad_type_str;
  string solver_type_str;
  string transport_correction_str;
                            
  // Read the type of the solver: SI or GMRES
  parameters_file>>solver_type_str;
  if (solver_type_str.compare("SI")==0 || solver_type_str.compare("si")==0)
    solver_type = si;
  else
  {
    if (solver_type_str.compare("BiCGSTAB")==0 ||
        solver_type_str.compare("BICGSTAB")==0 ||
        solver_type_str.compare("bicgstab")==0)
      solver_type = bicgstab;
    else
    {
      if (solver_type_str.compare("GMRES_CONDNUM")==0 ||
          solver_type_str.compare("gmres_condnum")==0)
        solver_type = gmres_condnum;
      else
        solver_type = gmres;
    }
  }

  // Read the tolerance on the solver
  parameters_file>>tolerance;
  
  // Read the maximum number of iterations
  parameters_file>>max_it;

  // Read the verbosity of the code
  parameters_file>>verbose;
  
  // Read if Fokker-Planck cross-section is used
  parameters_file>>fokker_planck_str;
  if (fokker_planck_str.compare("true")==0)
    fokker_planck = true;
  else
    fokker_planck = false;
  
  // Read if transport correction is used
  parameters_file>>transport_correction_str;
  if (transport_correction_str.compare("true")==0)
    transport_correction = true;
  else
    transport_correction = false;
  
  // Read if the optimal transport correction is used
  if (transport_correction==true)
  {
    parameters_file>>optimal_str;
    if (optimal_str.compare("true")==0)
      optimal = true;
    else
      optimal = false;
  }
  
  // Read if the angular multigrid preconditioning is used
  parameters_file>>multigrid_str;
  if (multigrid_str.compare("true")==0)
    multigrid = true;
  else 
    multigrid = false;
  
  // Read if MIP preconditioning is used
  parameters_file>>mip_str;
  if (mip_str.compare("true")==0)
    mip = true;
  else
    mip = false;

  // If MIP preconditioning is used, read the MIP solver type
  if (mip==true)
  {
    parameters_file>>mip_solver_type_str;
    if (mip_solver_type_str.compare("AGMG")==0 || 
        mip_solver_type_str.compare("agmg")==0)
      mip_solver_type = agmg;
    else
    {
      if (mip_solver_type_str.compare("CG_ML")==0 || 
          mip_solver_type_str.compare("cg_ml")==0)
      {
        string aggregation_type_str;
        mip_solver_type = cg_ml;
        // Read the aggregation type used
        parameters_file>>aggregation_type_str;
        if (aggregation_type_str.compare("Uncoupled")==0 ||
            aggregation_type_str.compare("uncoupled")==0)
          aggregation_type = uncoupled;
        else
        {
          if (aggregation_type_str.compare("MIS")==0 ||
              aggregation_type_str.compare("mis")==0)
            aggregation_type = mis;
          else
            aggregation_type = uncoupled_mis;
        }
      }
      else
      {
        if (mip_solver_type_str.compare("CG_SGS")==0 || 
            mip_solver_type_str.compare("cg_sgs")==0)
          mip_solver_type = cg_sgs;
        else
          mip_solver_type = cg_none;
      }
    }
  }
  
  // Read quadrature type: GLC (Gauss-Legendre-Chebyshev) or LS (Level
  // Symmetric)
  parameters_file>>quad_type_str;
  if (quad_type_str.compare("LS")==0)
    quad_type = ls;
  else
    quad_type = glc;
  
  // Read if the quadrature is a Galerkin quadrature
  parameters_file>>galerkin_str;
  if (galerkin_str.compare("true")==0)
    galerkin = true;
  else
    galerkin = false;
  
  // Read L_max
  parameters_file>>L_max;
  
  // Read order of the Sn order
  parameters_file>>sn;
  
  // Read the FEM type: BLD (Bilinear Discontinuous) of PWLD (PieceWise Linear
  // Discontinuous)
  parameters_file>>fe_type_str;
  if (fe_type_str.compare("BLD")==0)
    fe_type = bld;
  else
    fe_type = pwld;

  // Read the values of the source
  src.resize(n_src);
  for (unsigned int i=0; i<n_src; ++i)
    parameters_file>>src[i];
  // Read bottom incoming flux
  parameters_file>>inc_bottom;
  // Read right incoming flux
  parameters_file>>inc_right;
  // Read top incoming flux
  parameters_file>>inc_top;
  // Read left incoming
  parameters_file>>inc_left;
  // Read the total and the scattering cross section
  d_vector correction_vector;
  if (transport_correction==true && optimal==false)
    correction_vector.resize(n_mat,0.);
  sigma_t.resize(n_mat);
  unsigned int j_max(0);
  for (unsigned int jj=0; jj<=L_max; ++jj)
    j_max += jj+1;
  sigma_s.resize(n_mat,d_vector(j_max,0.));
  if (fokker_planck==true)
    alpha.resize(n_mat);
  for (unsigned int i=0; i<n_mat; ++i)
  {
    parameters_file>>sigma_t[i];
    if (fokker_planck==true)
      parameters_file>>alpha[i];
    else
    {                               
      for (unsigned int j=0; j<j_max; ++j)
        parameters_file>>sigma_s[i][j];
      if (transport_correction==true && optimal==false)
        parameters_file>>correction_vector[i];
    }
  }

  // Close the file
  parameters_file.close();

  Apply_parameters(n_mat,correction_vector);
}

void PARAMETERS::Apply_parameters(const unsigned int n_mat,
    d_vector const &correction_vector)
{
  // Create the Fokker-Planck cross-section
  if (fokker_planck==true)
    Build_fokker_planck_xs(n_mat);
  if (multigrid==false)
    n_levels = 1;
  else
    n_levels = ceil(log(double(sn))/log(2.));
  if (mip==true)
    ++n_levels;
  // Loop over the level of the angular multigrid
  sigma_t_lvl.resize(n_mat,d_vector(n_levels,0.));
  sigma_s_lvl.resize(n_mat,vector<d_vector>(n_levels));
  double L(L_max); 
  for (unsigned int lvl=0; lvl<n_levels; ++lvl)
  {
    for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
    {
      double L_2(ceil(L/2.));
      sigma_t_lvl[i_mat][lvl] = sigma_t[i_mat];
      sigma_s_lvl[i_mat][lvl] = d_vector(sigma_s[i_mat].begin(),
          sigma_s[i_mat].end());
      // Apply the transport correction
      if (transport_correction==true)
      {
        double correction(0);
        if (optimal==false)
        {
          if (lvl==0)
            correction = correction_vector[i_mat];
          else
            correction = sigma_s_lvl[i_mat][lvl-1][L_2]; 
        }
        Apply_transport_correction(i_mat,lvl,L_2,correction);
      }
    }
    L = ceil(L/2.);
    if (L==1.)
      L = 0.;
  }
}

void PARAMETERS::Build_fokker_planck_xs(const unsigned int n_mat)
{               
  for (unsigned int i_mat=0; i_mat<n_mat; ++i_mat)
  {
    unsigned int index(0);
    for (unsigned int l=0; l<L_max; ++l)
      for (unsigned int m=0; m<=l; ++m)
      {
        sigma_s[i_mat][index] = alpha[i_mat]/2.*(L_max*(L_max+1)-l*(l+1));
        ++index;
      }
  }
}

void PARAMETERS::Apply_transport_correction(unsigned int i_mat,unsigned int lvl,
    double L,double correction)
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
      correction = (sigma_s_lvl[i_mat][lvl][pos]+sigma_s_lvl[i_mat][lvl].back())/2.;
    }
  }
  
  sigma_t_lvl[i_mat][lvl] -= correction;
  d_vector::iterator sigma_s(sigma_s_lvl[i_mat][lvl].begin());
  d_vector::iterator sigma_s_end(sigma_s_lvl[i_mat][lvl].end());
  for (; sigma_s<sigma_s_end; ++sigma_s)
    *sigma_s -= correction;
}
