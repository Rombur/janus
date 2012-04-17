#ifndef _MIP_HH_
#define _MIP_HH

#include <cassert>
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
    MIP(unsigned int level,DOF_HANDLER* dof,PARAMETERS const* param,
        QUADRATURE const* quad,Epetra_Comm const* comm);

    ~MIP();

    /// Solve the MIP equation and add the solution to the flux moments.
    void Solve(Epetra_MultiVector &flux_moments);

    /// Delete the ML object before MPI_Finalize() is called.
    void Free_ml();

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

    /// Compute the right-hand side of the MIP equation.
    void Compute_rhs(Epetra_MultiVector const &x,Epetra_MultiVector &b);

    /// Build the CRS matrix of the MIP that will be used by CG and the
    /// algebraic multigrid method.
    void Build_lhs();

    /// Compute the penalty coefficient used by MIP.
    double Compute_penalty_coefficient(double D_m,double D_p,double h_m,double h_p,
        bool interior);

    /// Solve the system of equation using non-preconditioned Conjuguate
    /// Gradient.
    void Cg_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b);

    /// Solve the system of equation using Conjuguate Gradient preconditioned
    /// using Symmmetric-Gauss-Seidel.
    void Cg_sgs_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b);

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
    void Cg_ml_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b);

#ifdef AGMG
    /// Solve the system of equation using AGMG.
    void Agmg_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b);
#endif

    /// Add the solution to the 0th angular flux moment if the transport is
    /// solved using a Krylov method. Otherwise the correction is copied in
    /// flux_moments.
    void Project_solution(Epetra_MultiVector &flux_moments,
        Epetra_MultiVector const &x);
    
    /// Convert the Epetra left-hand side to fortran data types.
    void Convert_lhs_to_fortran(int* &ia,int* &ja,double* &a,unsigned int n_dof,
        unsigned int nnz);
    
    /// Convert the Epetra right-hand side to fortran data types.
    void Convert_rhs_to_fortran(Epetra_MultiVector const &b,double* f,
        const unsigned int n_dof) const;

    /// Convert the fortran rhs to Epetra rhs.
    void Convert_rhs_to_epetra(Epetra_MultiVector &b,double* f,unsigned int n_dof)    
      const;

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
