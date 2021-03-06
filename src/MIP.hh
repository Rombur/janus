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

#ifndef _MIP_HH_
#define _MIP_HH_

#include <cmath>
#include <iostream>
#include "gsl_math.h"
#include "AztecOO.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Ifpack_PointRelaxation.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "CELL.hh"
#include "EDGE.hh"
#include "ERROR_ESTIMATOR.hh"
#include "EXCEPTION.hh"
#include "DOF_HANDLER.hh"
#include "FINITE_ELEMENT.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"

using namespace std;

#ifdef AGMG
extern "C"
{
  /**
   * Wrapper for the Fortran function dagmg which solve the problem using an
   * algebraic multigrid method.
   * @param n order of the linear system
   * @param a nonzero matrix entries (numerical values), stored row-wise
   * @param ja corresponding column entries
   * @param ia indicate where every row starts
   * @param f right-hand side
   * @param x initial guess (optional) and output
   * @param ijob 0 = setup+solve+memory release, no initial guess
   * @param iprint unit number for output messages (6 is the screen)
   * @param nrest 1 = flexible CG
   * @param iter maximum number of iteration on input and effective number of
   * iteration on output
   * @param tol = tolerance on the relative residual norm (||r||_2/||b||_2)
   * @todo use dagmgpar_
   * @todo use setup and then solve
   */
  void dagmg_(int* n,double a[],int ja[],int ia[],double f[],double x[],int* ijob,
      int* iprint,int* nrest,int* iter,double* tol);
}
#endif

/**
 * This class builds and solves the MIP DSA. The matrix is built using the CRS
 * format and then solved using CG preconditioned by symmetric-Gauss-Seidel,
 * ML or AGMG.
 */

class MIP
{
  public :
    MIP(Epetra_Comm const* comm,PARAMETERS const* param,DOF_HANDLER* dof,
        QUADRATURE const* quad=NULL,unsigned int level=0);

    ~MIP();

    /// Delete the ML object before MPI_Finalize() is called.
    void Free_ml();

    /// Set the value of group and if necessary delete A or a and the
    /// preconditioner for CG.
    void Set_group(unsigned int g);

    /// Solve the MIP equations and add the solution to the flux moments.
    void Solve(Epetra_MultiVector &flux_moments,
        Epetra_MultiVector* initial_guess=NULL);

    /// Solve the multigroup MIP equations when MIP is used alone.
    void Solve_diffusion(const unsigned int n_groups,Epetra_MultiVector &flux_moments,
        Epetra_MultiVector const &group_flux,Epetra_MultiVector* initial_guess=NULL);

    /// Return the time needed to build the matrix.
    double Get_building_mip_time() const;

    /// Return the time need to initialize the preconditioner for CG
    double Get_init_prec_mip_time() const;

    /// Return the time needed to solve the matrix. When using AGMG, the
    /// conversion time between data types is not counted.
    double Get_solve_mip_time() const;

  private :
    /// Compute the number of non-zero elements for each row of #A.
    void Compute_n_entries_per_row(int* n);

    /// Compute the right-hand side of the MIP equations.
    void Compute_rhs(Epetra_MultiVector const &x,Epetra_MultiVector &b);

    /// Compute the right-hand side of the multigroup MIP equations when MIP
    /// is used.
    void Compute_diffusion_rhs(const unsigned int ngroups,Epetra_MultiVector const &x,
        Epetra_MultiVector &b);

    /// Build the CRS matrix of the MIP that will be used by CG and the
    /// algebraic multigrid method.
    void Build_lhs();

    /// Compute the penalty coefficient used by MIP.
    double Compute_penalty_coefficient(double D_m,double D_p,double h_m,double h_p,
        bool interior);

    /// Solve the system of equation using non-preconditioned Conjuguate
    /// Gradient.
    void Cg_solve(const bool dsa,Epetra_MultiVector &flux_moments,
        Epetra_MultiVector &b,Epetra_MultiVector* initial_guess);

    /// Solve the system of equation using Conjuguate Gradient preconditioned
    /// using SSOR.
    void Cg_ssor_solve(const bool dsa,Epetra_MultiVector &flux_moments,
        Epetra_MultiVector &b,Epetra_MultiVector* initial_guess);

    /**
     * Solve the system of equation using Conjuguate Gradient preconditioned
     * using ML. 
     * The default paramaters are:
     *  - "max levels" = 10
     *  - "prec type" = "MGV"
     *  - "increasing or decreasing" = "increasing"
     *  - "aggregation: type" = "Uncoupled-MIS"
     *  - "aggregation: damping factor" = 1.333
     *  - "eigen-analysis: type" = "cg"
     *  - "eigen-analysis: iterations" = 10
     *  - "smoother: sweeps" = 2
     *  - "smoother: damping factor" = 1.0
     *  - "smoother: pre or post" = "both"
     *  - "smoother: type" = "symmetric Gauss-Seidel"
     *  - "coarse: type" = "Amesos-KLU"
     *  - "coarse: max size" = 128
     */
    void Cg_ml_solve(const bool dsa,Epetra_MultiVector &flux_moments,
        Epetra_MultiVector &b,Epetra_MultiVector* initial_guess);

#ifdef AGMG
    /// Solve the system of equation using AGMG.
    void Agmg_solve(const bool dsa,Epetra_MultiVector &flux_moments,
        Epetra_MultiVector &b,Epetra_MultiVector* initial_guess);
#endif

    /// Add the solution to the 0th angular flux moment if the transport is
    /// solved using a Krylov method. Otherwise the correction is copied in
    /// flux_moments.
    void Project_solution(Epetra_MultiVector &flux_moments,
        Epetra_MultiVector const &x);
    
    /// Convert the Epetra left-hand side to fortran data types.
    void Convert_lhs_to_fortran(int* &ia,int* &ja,double* &a,unsigned int n_dof,
        unsigned int nnz);
    
    /// Convert the Epetra MultiVector to a fortran data type.
    void Convert_to_fortran(Epetra_MultiVector const &b,double* f,
        const unsigned int n_dof) const;

    /// Convert the fortran data type to Epetra MultiVector.
    void Convert_to_epetra(Epetra_MultiVector &b,double* f,unsigned int n_dof)    
      const;

    /// Current group.
    unsigned int group;
    /// Level of the cross sections to use.
    unsigned int lvl;
    /// Indicate where every row starts.
    int* ia;
    /// Column entries.
    int* ja;
    /// Nonzero matrix entries (numerical values), stored row-wise.
    double* a;
    /// Epetra communicator.
    Epetra_Comm const* comm;
    /// Pointer to the dof_handler object. (DOF_HANDLER is not const because
    /// we need to return an iterator on the cells).
    DOF_HANDLER* dof_handler;
    /// Pointer to the parameters object.
    PARAMETERS const* parameters;
    /// Pointer to the quadrature object.
    QUADRATURE const* quad;
    /// Epetra map for the left-hand side matrix.
    Epetra_Map* mip_map;
    /// Left-hand side CRS matrix of MIP. The matrix is built the first time
    /// that #Solve is called and then is reused.
    Epetra_CrsMatrix* A;
    /// Timer for building the matrix.
    Teuchos::Time* building_timer;
    /// Timer for initializing the CG preconditioner.
    Teuchos::Time* init_prec_timer;
    /// Timer for solving MIP.
    Teuchos::Time* solve_timer;
    /// SGS preconditioner.
    Ifpack_PointRelaxation* sgs_prec;
    /// Multi-level preconditioner.
    ML_Epetra::MultiLevelPreconditioner* ml_prec;
};

inline double MIP::Get_building_mip_time() const
{
  return building_timer->totalElapsedTime();
}

inline double MIP::Get_init_prec_mip_time() const
{
  return init_prec_timer->totalElapsedTime();
}

inline double MIP::Get_solve_mip_time() const
{
  return solve_timer->totalElapsedTime();
}

#endif 
