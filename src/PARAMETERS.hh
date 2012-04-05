#ifndef _PARAMETERS_HH_
#define _PARAMETERS_HH_

#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <iostream>

using namespace std;

typedef vector<double> d_vector;

/// Enum on the aggregation type used by ML: uncoupled, mis or uncoupled-mis
enum AGGREGATION_TYPE{uncoupled,mis,uncoupled_mis};

/// Enum on the FE type: bld (BiLinear Discontinuous) and pwld (PieceWise
/// Linear Discontinuous).
enum FE_TYPE{bld,pwld};

/// Enum on the MIP solver type: cg_none (CG without preconditioning), cg_sgs
/// (CG preconditioned with Symmetric-Gauss-Seidel), cg_ml (CG preconditioned
/// with ML) and agmg (AGMG).
enum MIP_SOLVER_TYPE{cg_none,cg_sgs,cg_ml,agmg};

/// Enum on the quadrature type: glc (Gauss-Legendre-Chebyshev) or ls (Level
/// Symmetric).
enum QUAD_TYPE{glc,ls};

/// Enum on the solver type: si (Source Iteration), gmres, gmres_condnum
/// (GMRES with estimation of the condition number) or bicgstab (BiCGSTAB).
enum SOLVER_TYPE{si,gmres,gmres_condnum,bicgstab};

/**
 * Read and store all the parameters. Parameters are read in the following
 * order: solver type, tolerance, maximum number of iterations, Fokker-Planck
 * flag, transport correction flag, optimal transport correction flag,
 * multigrid flag, mip flag, quadrature type, Galerkin flag, L_max, Sn order,
 * fe type, intensities of the source, bottom incoming flux, right incoming
 * flux, top incoming flux, left incoming flux, sigma_t, sigma_s.
 */
class PARAMETERS
{
  public :
    PARAMETERS(string* parameters_inputfile);

    /// Read the parameters input file.
    void Read_parameters(unsigned int n_src,unsigned int n_mat);

    /// Return the flag on Galerkin quadrature.
    bool Get_galerkin() const;

    /// Return true if angular multigrid is used.
    bool Get_multigrid() const;

    /// Return true if MIP is used.
    bool Get_mip() const;

    /// Return the maximum number of iterations.
    unsigned int Get_max_it() const;

    /// Return #L_max.
    unsigned int Get_L_max() const;

    /// Return Sn order.
    unsigned int Get_sn_order() const;

    /// Return the number of "angular levels" (DSA counts as an angular level).
    unsigned int Get_n_levels() const;

    /// Return the level of verbosity of the code.
    unsigned int Get_verbose() const;

    /// Return the relative tolerance.
    double Get_tolerance() const;

    /// Return the incoming flux of the bottom boundary.
    double Get_inc_bottom() const;

    /// Return the incoming flux of the right boundary.
    double Get_inc_right() const;

    /// Return the incoming flux of the top boundary.
    double Get_inc_top() const;

    /// Return the incoming flux of the left boundary.
    double Get_inc_left() const;

    /// Return the intensity of the source i.
    double Get_src(unsigned int i) const;

    /// Return the type iof aggregation used by ML (uncoupled, mis or
    /// uncoupled-mis).
    AGGREGATION_TYPE Get_aggregation_type() const;

    /// Return the type of finite elements used (BLD or PWLD).
    FE_TYPE Get_fe_type() const;

    /// Return the type of quadrature used (GLC or LS). 
    QUAD_TYPE Get_quad_type() const;

    /// Return the type of solver used (SI, BiCGSTAB, GMRES or GMRES with
    /// estimation of the condition number).
    SOLVER_TYPE Get_solver_type() const;

    /// Return the type of MIP solver used (CG without preconditioning, CG
    /// with SGS preconditioning, CG with ML preconditioning or AGMG).
    MIP_SOLVER_TYPE Get_mip_solver_type() const;

    /// Return sigma_t for the material i.
    d_vector Get_sigma_t(unsigned int i) const;

    /// Return sigma_s for the material i.
    vector<d_vector> Get_sigma_s(unsigned int i) const;

  private :
    /// Apply all the parameters.
    void Apply_parameters(const unsigned int n_mat,d_vector const &correction_vector);

    /// Build the Fokker-Planck cross sections.
    void Build_fokker_planck_xs(const unsigned int n_mat);

    /// Apply the "standard" extended transport correction or the optimal
    /// transport corrections.
    void Apply_transport_correction(unsigned int i_mat,unsigned int lvl,double L,
        double correction);

    /// If flag is true, Fokker-Planck cross-section \f$\left(\frac{\alpha}{2} (L(L+1)-l(l+1))\right)\f$ is used.
    bool fokker_planck;
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
    /// Maximum number of iterations.
    unsigned int max_it;
    /// L_max.
    unsigned int L_max;
    /// Sn order.
    unsigned int sn;
    /// Number of levels of the angular multigrid.
    unsigned int n_levels;
    /// Level of verbosity of the code.
    unsigned int verbose;
    /// Tolerance on the solver.
    double tolerance;
    /// Bottom incoming flux.
    double inc_bottom;
    /// Right incoming flux.
    double inc_right;
    /// Top incoming flux.
    double inc_top;
    /// Left incoming flux.
    double inc_left;
    /// Type of aggregation used by ML: uncoupled, mis or uncoupld-mis.
    AGGREGATION_TYPE aggregation_type;
    /// Type of finite elements: BLD or PWLD.
    FE_TYPE fe_type;
    /// Type of quadrature: GLC or LS.
    QUAD_TYPE quad_type;
    /// Type of MIP solver type: CG without preconditioning, CG preconditioned with 
    /// Symmetric-Gauss-Seidel, CG preconditioned with ML and AGMG.
    MIP_SOLVER_TYPE mip_solver_type;
    /// Type of solver: SI, BiCGSTAB or GMRES.
    SOLVER_TYPE solver_type;
    /// Pointer to the name of the inpuit file containing the parameters of
    /// the problem.
    string* parameters_filename;
    /// Values of the source.
    d_vector src;
    /// Total cross section.
    d_vector sigma_t;
    /// Scattering cross section.
    vector<d_vector> sigma_s;
    /// \f$\alpha\f$ parameter of the Fokker-Planck cross sections.
    d_vector alpha;
    /// Total cross section for the different level of the angular multigrid.
    vector<d_vector> sigma_t_lvl;
    /// Scattering cross section for the different level of the angular
    /// multigrid.
    vector<vector<d_vector> > sigma_s_lvl;
};

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

inline unsigned int PARAMETERS::Get_max_it() const
{
  return max_it;
}

inline unsigned int PARAMETERS::Get_L_max() const
{
  return L_max;
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

inline double PARAMETERS::Get_tolerance() const
{
  return tolerance;
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

inline double PARAMETERS::Get_src(unsigned int i) const
{
  return src[i];
}
inline AGGREGATION_TYPE PARAMETERS::Get_aggregation_type() const
{
  return aggregation_type;
}

inline FE_TYPE PARAMETERS::Get_fe_type() const
{
  return fe_type;
}

inline QUAD_TYPE PARAMETERS::Get_quad_type() const
{
  return quad_type;
}

inline MIP_SOLVER_TYPE PARAMETERS::Get_mip_solver_type() const
{
  return mip_solver_type;
}

inline SOLVER_TYPE PARAMETERS::Get_solver_type() const
{
  return solver_type;
}

inline d_vector PARAMETERS::Get_sigma_t(unsigned int i) const
{
  return sigma_t_lvl[i];
}

inline vector<d_vector> PARAMETERS::Get_sigma_s(unsigned int i) const
{
  return sigma_s_lvl[i];
}

#endif
