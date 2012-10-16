/*
Copyright (c) 2012, Bruno Turcksin.

This file is part of Janus.

Janu is free software: you can redistribute it and/or modify
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

#include "PARAMETERS.hh"

PARAMETERS::PARAMETERS(string* parameters_inputfile) :
  parameters_filename(parameters_inputfile)
{}

void PARAMETERS::Read_parameters(const unsigned int n_src)
{
  // Open the file to read it
  ifstream parameters_file(parameters_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(parameters_file.good(),string ("unable to open the file "+
        *parameters_filename + " containing the parameters."));

  string weight_sum_str;
  string bottom_bc_type_str;
  string left_bc_type_str;
  string right_bc_type_str;
  string top_bc_type_str;
  string fe_type_str;
  string xs_type_str;
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

  // Read the tolerance for the inner solver
  parameters_file>>inner_tolerance;

  // Read the tolerance for the outer solvers
  parameters_file>>group_tolerance;
  
  // Read the maximum number of inner iterations
  parameters_file>>max_inner_it;
  
  // Read the maximum number of supergroup iterations
  parameters_file>>max_supergroup_it;
  
  // Read the maximum number of group iterations
  parameters_file>>max_group_it;

  // Read the sum of the weight for the quadrature
  parameters_file>>weight_sum_str;
  if (weight_sum_str.compare("2_PI")==0 || weight_sum_str.compare("2_pi")==0)
    weight_sum = 2.*M_PI;
  else
  {
    if (weight_sum_str.compare("4_PI")==0 || weight_sum_str.compare("4_pi")==0)
      weight_sum = 4.*M_PI;
    else
      weight_sum = 1.0;
  }
  
  // Read the type of cross section file
  parameters_file>>xs_type_str;
  if (xs_type_str.compare("fp")==0)
    xs_type = fp;
  else
  {
    if (xs_type_str.compare("regular")==0)
      xs_type = regular;
    else
    {
      if (xs_type_str.compare("regular_exs")==0)
        xs_type = regular_exs;
      else
      {
        Check(xs_type_str.compare("cepxs")==0,string ("Unknown cross section type."));
        xs_type = cepxs;
      }
    }
    string permutation_type_str;
    parameters_file>>permutation_type_str;
    if (permutation_type_str.compare("none")==0)
      permutation_type = none;
    else
    {
      if (permutation_type_str.compare("logarithmic")==0)
        permutation_type = logarithmic;
      else
        permutation_type = linear;
    }
  }

  // Read the verbosity of the code
  parameters_file>>verbose;
  
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
        if (mip_solver_type_str.compare("CG_SSOR")==0 || 
            mip_solver_type_str.compare("cg_ssor")==0)
        {

          mip_solver_type = cg_ssor;
          parameters_file>>damping_factor;
        }
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
  {
    Check(quad_type_str.compare("GLC")==0,string ("Unknown quadrature type."));
    quad_type = glc;
  }
  
  // Read if the quadrature is a Galerkin quadrature
  parameters_file>>galerkin_str;
  if (galerkin_str.compare("true")==0)
    galerkin = true;
  else
    galerkin = false;
  
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
  unsigned int n_groups(0);
  parameters_file>>n_groups;
  src.resize(n_src,d_vector (n_groups));
  for (unsigned int i=0; i<n_src; ++i)
    for (unsigned int g=0; g<n_groups; ++g)
      parameters_file>>src[i][g];
  
  // Read the type of boundary condition on the bottom side
  parameters_file>>bottom_bc_type_str;
  if (bottom_bc_type_str.compare("vacuum")==0)
    bottom_bc_type = vacuum;
  else
  {
    if (bottom_bc_type_str.compare("reflective")==0)
      bottom_bc_type = reflective;
    else
    {
      if (bottom_bc_type_str.compare("most_normal")==0)
        bottom_bc_type = most_normal;
      else
      {
        Check(bottom_bc_type_str.compare("isotropic")==0,
            string ("Unknown boundary condition type on the bottom boundary."));
        bottom_bc_type = isotropic;
      }
      // Read bottom incoming flux
      parameters_file>>inc_bottom;
    }
  }
  // Read the type of boundary condition on the right side
  parameters_file>>right_bc_type_str;
  if (right_bc_type_str.compare("vacuum")==0)
    right_bc_type = vacuum;
  else
  {
    if (right_bc_type_str.compare("reflective")==0)
      right_bc_type = reflective;
    else
    {
      if (right_bc_type_str.compare("most_normal")==0)
        right_bc_type = most_normal;
      else
      {
        Check(right_bc_type_str.compare("isotropic")==0,
            string ("Unknown boundary condition type on the right boundary."));
        right_bc_type = isotropic;
      }
      // Read right incoming flux
      parameters_file>>inc_right;
    }
  }
  // Read the type of boundary condition on the top side
  parameters_file>>top_bc_type_str;
  if (top_bc_type_str.compare("vacuum")==0)
    top_bc_type = vacuum;
  else
  {
    if (top_bc_type_str.compare("reflective")==0)
      top_bc_type = reflective;
    else
    {
      if (top_bc_type_str.compare("most_normal")==0)
        top_bc_type = most_normal;
      else
      {
        Check(top_bc_type_str.compare("isotropic")==0,
            string ("Unknown boundary condition type on the top boundary."));
        top_bc_type = isotropic;
      }
      // Read right incoming flux
      parameters_file>>inc_top;
    }
  }
  // Read the type of boundary condition on the left side
  parameters_file>>left_bc_type_str;
  if (left_bc_type_str.compare("vacuum")==0)
    left_bc_type = vacuum;
  else
  {
    if (left_bc_type_str.compare("reflective")==0)
      left_bc_type = reflective;
    else
    {
      if (left_bc_type_str.compare("most_normal")==0)
        left_bc_type = most_normal;
      else
      {
        Check(left_bc_type_str.compare("isotropic")==0,
            string ("Unknown boundary condition type on the left boundary."));
          left_bc_type = isotropic;
      }
      // Read right incoming flux
      parameters_file>>inc_left;
    }
  }

  // Close the file
  parameters_file.close();

  // Compute the number of level when the angular multigrid is used
  if (multigrid==false)
    n_levels = 1;
  else
    n_levels = ceil(log(double(sn))/log(2.));
}
