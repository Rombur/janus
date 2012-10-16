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

#ifndef _PARAMETERS_HH_
#define _PARAMETERS_HH_

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "gsl_math.h"
#include "EXCEPTION.hh"

using namespace std;

typedef vector<double> d_vector;

/// Enum on the aggregation type used by ML: uncoupled, mis or uncoupled-mis
enum AGGREGATION_TYPE{uncoupled,mis,uncoupled_mis};

/// Enum on the different types of boundary conditions.
enum BC_TYPE{vacuum,isotropic,most_normal,reflective};

/// Enum on the FE type: bld (BiLinear Discontinuous) and pwld (PieceWise
/// Linear Discontinuous).
enum FE_TYPE{bld,pwld};

/// Enum on the MIP solver type: cg_none (CG without preconditioning), cg_ssor
/// (CG preconditioned with SSOR), cg_ml (CG preconditioned with ML), and agmg (AGMG).
enum MIP_SOLVER_TYPE{cg_none,cg_ssor,cg_ml,agmg};

/// Enum on the permutation type used with CEPXS cross sections: none, linear
/// or logarithmic.
enum PERMUTATION_TYPE{none,linear,logarithmic};

/// Enum on the quadrature type: glc (Gauss-Legendre-Chebyshev) or ls (Level
/// Symmetric).
enum QUAD_TYPE{glc,ls};

/// Enum on the solver type: si (Source Iteration), gmres, gmres_condnum
/// (GMRES with estimation of the condition number) or bicgstab (BiCGSTAB).
enum SOLVER_TYPE{si,gmres,gmres_condnum,bicgstab};

/// Enum on the type of cross section file: fp (Fokker-Planck), regular,
/// regular_exs (regular with energy deposition cross sections), cepxs
/// (CEPXS).
enum XS_TYPE{fp,regular,regular_exs,cepxs};

/**
 * Read and store all the parameters. Parameters are read in the following
 * order: solver type, tolerance, maximum number of iterations, Fokker-Planck
 * flag, transport correction flag, optimal transport correction flag,
 * multigrid flag, mip flag, quadrature type, Galerkin flag, Sn order,
 * fe type, intensities of the source, bottom incoming flux, right incoming
 * flux, top incoming flux, left incoming flux.
 */
class PARAMETERS
{
  public :
    PARAMETERS(string* parameters_inputfile);

    /// Read the parameters input file.
    void Read_parameters(const unsigned int n_src);

    /// Return the flag on the transport correction.
    bool Get_transport_correction() const;

    /// Return the flag on the optimal transport correction.
    bool Get_optimal_tc() const;

    /// Return the flag on Galerkin quadrature.
    bool Get_galerkin() const;

    /// Return true if angular multigrid is used.
    bool Get_multigrid() const;

    /// Return true if MIP is used.
    bool Get_mip() const;

    /// Return the maximum number of inner iterations.
    unsigned int Get_max_inner_it() const;

    /// Return the maximum number of iterations on each supergroup.
    unsigned int Get_max_supergroup_it() const;

    /// Return the maximum number of iterations on all the groups.
    unsigned int Get_max_group_it() const;

    /// Return Sn order.
    unsigned int Get_sn_order() const;

    /// Return the number of "angular levels" (DSA counts as an angular
    /// level).
    inline unsigned int Get_n_levels() const;

    /// Return the level of verbosity of the code.
    unsigned int Get_verbose() const;

    /// Return the relative tolerance for the inner solver.
    double Get_inner_tolerance() const;

    /// Return the relative tolerance for the group and supergroup solver.
    double Get_group_tolerance() const;

    /// Return the incoming flux of the bottom boundary.
    double Get_inc_bottom() const;

    /// Return the incoming flux of the right boundary.
    double Get_inc_right() const;

    /// Return the incoming flux of the top boundary.
    double Get_inc_top() const;

    /// Return the incoming flux of the left boundary.
    double Get_inc_left() const;

    /// Return the intensity of the source i.
    d_vector Get_src(unsigned int i) const;

    /// Return the sum of the weights (1, 2pi or 4pi)
    double Get_weight_sum() const;

    /// Return the damping factor used by SSOR.
    double Get_damping_factor() const;

    /// Return the type of aggregation used by ML (uncoupled, mis or
    /// uncoupled-mis).
    AGGREGATION_TYPE Get_aggregation_type() const;

    /// Return the boundary condition of the bottom side of the domain
    /// (vacuum, isotropic, most_normal or reflective).
    BC_TYPE Get_bottom_bc_type() const;

    /// Return the boundary condition of the right side of the domain
    /// (vacuum, isotropic, most_normal or reflective).
    BC_TYPE Get_right_bc_type() const;

    /// Return the boundary condition of the top side of the domain
    /// (vacuum, isotropic, most_normal or reflective).
    BC_TYPE Get_top_bc_type() const;

    /// Return the boundary condition of the left side of the domain
    /// (vacuum, isotropic, most_normal or reflective).
    BC_TYPE Get_left_bc_type() const;

    /// Return the type of finite elements used (BLD or PWLD).
    FE_TYPE Get_fe_type() const;

    /// Return the type of MIP solver used (CG without preconditioning, CG
    /// with SGS preconditioning, CG with ML preconditioning or AGMG).
    MIP_SOLVER_TYPE Get_mip_solver_type() const;

    /// Return the type of permutation used when CEPXS cross sections are
    /// used: none, linear or logarithmic.
    PERMUTATION_TYPE Get_permutation_type() const;

    /// Return the type of quadrature used (GLC or LS). 
    QUAD_TYPE Get_quad_type() const;

    /// Return the type of solver used (SI, BiCGSTAB, GMRES or GMRES with
    /// estimation of the condition number).
    SOLVER_TYPE Get_solver_type() const;

    /// Return the type of cross section file used (Fokker-Planck, regular,
    /// regular with energy deposition cross sections, CEPXS).
    XS_TYPE Get_xs_type() const;

  private :
    /// If flag is true, the transport correction is used.
    bool transport_correction;
    /// If flag is true and that #transport_correction is true, the optimal
    /// transport correction is used. Otherwise, the "standard" extended cross
    /// section is used.
    bool optimal;
    /// If flag is true, the angular multigrid preconditioning is used. GMRES
    /// has to be used.
    bool multigrid;
    /// If flag is true, MIP DSA is used.
    bool mip;
    /// If flas is true, a Galerkin quadrature is used.
    bool galerkin;
    /// Maximum number of inner iterations.
    unsigned int max_inner_it;
    /// Maximum number of iterations for each supergroup.
    unsigned int max_supergroup_it;
    /// Maximum number of iterations over all the groups.
    unsigned int max_group_it;
    /// Sn order.
    unsigned int sn;
    /// Number of levels of the angular multigrid.
    unsigned int n_levels;
    /// Level of verbosity of the code.
    unsigned int verbose;
    /// Tolerance on the inner solver (SI or Krylov solver).
    double inner_tolerance;
    /// Tolerance on the outer solver for the groups and supergroups.
    double group_tolerance;
    /// Bottom incoming flux.
    double inc_bottom;
    /// Right incoming flux.
    double inc_right;
    /// Top incoming flux.
    double inc_top;
    /// Left incoming flux.
    double inc_left;
    /// Sum of the weights (1, 2pi or 4pi)
    double weight_sum;
    /// Damping factor used by SSOR.
    double damping_factor;
    /// Type of aggregation used by ML: uncoupled, mis or uncoupld-mis.
    AGGREGATION_TYPE aggregation_type;
    /// Type of boundary condition on the bottom side.
    BC_TYPE bottom_bc_type;
    /// Type of boundary condition on the right side.
    BC_TYPE right_bc_type;
    /// Type of boundary condition on the top side.
    BC_TYPE top_bc_type;
    /// Type of boundary condition on the left side.
    BC_TYPE left_bc_type;
    /// Type of finite elements: BLD or PWLD.
    FE_TYPE fe_type;
    /// Type of quadrature: GLC or LS.
    QUAD_TYPE quad_type;
    /// Type of MIP solver type: CG without preconditioning, CG preconditioned with 
    /// Symmetric-Gauss-Seidel, CG preconditioned with ML and AGMG.
    MIP_SOLVER_TYPE mip_solver_type;
    /// Type of permutation used when CEPXS cross sections are used.
    PERMUTATION_TYPE permutation_type;
    /// Type of solver: SI, BiCGSTAB or GMRES.
    SOLVER_TYPE solver_type;
    /// Type of cross section file: Fokker-Planck, regular, regular with
    /// energy deposition cross section, CEPXS.
    XS_TYPE xs_type;
    /// Pointer to the name of the input file containing the parameters of
    /// the problem.
    string* parameters_filename;
    /// Values of the source.
    vector<d_vector> src;
};

inline bool PARAMETERS::Get_transport_correction() const
{
  return transport_correction;
}

inline bool PARAMETERS::Get_optimal_tc() const
{
  return optimal;
}

inline bool PARAMETERS::Get_galerkin() const
{
  return galerkin;
}

inline bool PARAMETERS::Get_multigrid() const
{
  return multigrid;
}

inline bool PARAMETERS::Get_mip() const
{
  return mip;
}

inline unsigned int PARAMETERS::Get_max_inner_it() const
{
  return max_inner_it;
}

inline unsigned int PARAMETERS::Get_max_supergroup_it() const
{
  return max_supergroup_it;
}

inline unsigned int PARAMETERS::Get_max_group_it() const
{
  return max_group_it;
}

inline unsigned int PARAMETERS::Get_sn_order() const
{
  return sn;
}

inline unsigned int PARAMETERS::Get_n_levels() const
{
  return n_levels;
}

inline unsigned int PARAMETERS::Get_verbose() const
{
  return verbose;
}

inline double PARAMETERS::Get_inner_tolerance() const
{
  return inner_tolerance;
}

inline double PARAMETERS::Get_group_tolerance() const
{
  return group_tolerance;
}

inline double PARAMETERS::Get_inc_bottom() const
{
  return inc_bottom;
}

inline double PARAMETERS::Get_inc_right() const
{
  return inc_right;
}

inline double PARAMETERS::Get_inc_top() const
{
  return inc_top;
}

inline double PARAMETERS::Get_inc_left() const
{
  return inc_left;
}

inline d_vector PARAMETERS::Get_src(unsigned int i) const
{
  return src[i];
}

inline double PARAMETERS::Get_weight_sum() const
{
  return weight_sum;
}

inline double PARAMETERS::Get_damping_factor() const
{
  return damping_factor;
}

inline AGGREGATION_TYPE PARAMETERS::Get_aggregation_type() const
{
  return aggregation_type;
}

inline BC_TYPE PARAMETERS::Get_bottom_bc_type() const
{
  return bottom_bc_type;
}

inline BC_TYPE PARAMETERS::Get_right_bc_type() const
{
  return right_bc_type;
}

inline BC_TYPE PARAMETERS::Get_top_bc_type() const
{
  return top_bc_type;
}

inline BC_TYPE PARAMETERS::Get_left_bc_type() const
{
  return left_bc_type;
}

inline FE_TYPE PARAMETERS::Get_fe_type() const
{
  return fe_type;
}

inline MIP_SOLVER_TYPE PARAMETERS::Get_mip_solver_type() const
{
  return mip_solver_type;
}

inline PERMUTATION_TYPE PARAMETERS::Get_permutation_type() const
{
  return permutation_type;
}

inline QUAD_TYPE PARAMETERS::Get_quad_type() const
{
  return quad_type;
}

inline SOLVER_TYPE PARAMETERS::Get_solver_type() const
{
  return solver_type;
}

inline XS_TYPE PARAMETERS::Get_xs_type() const
{
  return xs_type;
}

#endif
