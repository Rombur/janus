#ifndef _FINITE_ELEMENT_HH_
#define _FINITE_ELEMENT_HH_

#include <vector>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "CELL.hh"

using namespace std;

typedef vector<double> d_vector;

/**
 * Define the API for all the finite element discretizations.
 */ 

class FINITE_ELEMENT
{
  public :
    FINITE_ELEMENT(d_vector const &cell_x,d_vector const &cell_y) : 
      x(cell_x),y(cell_y){};

    /// Purely virtual function to build the 1-dimensional finite elements.
    virtual void Build_fe_1d() = 0;

    /// Purely virtual function to build the 2-dimensional finite elements.
    virtual void Build_fe_2d() = 0;

    /// Purely virtual function to build the upwinding matrices.
    virtual void Build_upwind_matrices(CELL &cell,vector<CELL> const &mesh) = 0;

    /// Return the number of degrees of freedom for the cell.
    unsigned int Get_dof_per_cell() const;

  protected :
    /// Number of degrees of freedom for the cell.
    unsigned int dof_per_cell;
    /// Abscissae of the cell associated to the element.
    d_vector x;
    /// Ordinates of the cell associated to the element.
    d_vector y;
    /// Downwind matrices \f$\int_E b_i\ b_j\ dr \f$ where \f$ b_i \f$ and
    /// \f$ b_j \f$ are defined on the same cell. 
    vector<Teuchos::SerialDenseMatrix<int,double> > downwind;
    /// Upwind matrices \f$\int_{E_c} b_i\ b_j\ dr \f$ where \f$ b_i \f$ and
    /// \f$ b_j \f$ are defined on different cell.
    vector<Teuchos::SerialDenseMatrix<int,double> > upwind;
    /// Compute \f$\int_E \partial_n b_i\ b_j\ dr \f$ \f$ where b_i$ and $ b_j\f$
    /// are defined on the same cell.
    vector<vector<Teuchos::SerialDenseMatrix<int,double> > > edge_deln_matrix;
    /// Compute \f$\int_{E_c} \partial_n b_i\ b_j\ dr \f$ \f$ where b_i$ and $ b_j\f$
    /// are defined on different cells.
    vector<vector<Teuchos::SerialDenseMatrix<int,double> > > coupling_edge_deln_matrix;
    /// Compute the mass matrix \f$\int_D b_i\ b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> mass_matrix;
    /// Compute the x component of the gradient matrix \f$\int_D b_i\ \partial_x b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> x_grad_matrix;
    /// Compute the x component of the gradient matrix \f$\int_D b_i\ \partial_y b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> y_grad_matrix;
    /// Compute the stiffness matrix \f$\int_D \boldsymbol{\nabla} b_i\ \boldsymbol{\nabla} b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> stiffness_matrix;
};

inline unsigned int FINITE_ELEMENT::Get_dof_per_cell() const
{
  return dof_per_cell;
}

#endif
