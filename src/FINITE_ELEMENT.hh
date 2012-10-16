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

#ifndef _FINITE_ELEMENT_HH_
#define _FINITE_ELEMENT_HH_

#include <vector>
#include "Teuchos_SerialDenseMatrix.hpp"
//#include "CELL.hh"

using namespace std;

typedef vector<double> d_vector;

class CELL;
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

    /// Purely virtual function to build the upwinding matrices.
    virtual void Build_upwind_matrices(CELL* cell,vector<CELL*> const &mesh) = 0;

    /// Purely virtual function to build the 2-dimensional finite elements.
    virtual void Build_fe_2d() = 0;

    /// Return the number of degrees of freedom for the cell.
    unsigned int Get_dof_per_cell() const;

    /// Return a pointer to the #downwind matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_downwind_matrix(
        unsigned int i) const;

    /// Return a pointer the #upwind matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_upwind_matrix(
        unsigned int i) const;

    /// Return a pointer to #edge_deln_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_edge_deln_matrix(
        unsigned int i,unsigned int j) const;

    /// Return a pointer to #coupling_edge_deln_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_coupling_edge_deln_matrix(
        unsigned int i,unsigned int j) const;

    /// Return a pointer to #mass_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_mass_matrix() const;

    /// Return a pointer to #x_grad_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_x_grad_matrix() const;

    /// Return a pointer to #y_grad_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_y_grad_matrix() const;

    /// Return a pointer to #x_grad_edge_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_x_grad_matrix(
        unsigned int i) const;

    /// Return a pointer to #y_grad_edge_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_y_grad_matrix(
        unsigned int i) const;

    /// Return a pointer to #stiffness_matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_stiffness_matrix() const;

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
    /// Mass matrix \f$\int_D b_i\ b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> mass_matrix;
    /// x component of the gradient matrix \f$\int_D b_i\ \partial_x b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> x_grad_matrix;
    /// y component of the gradient matrix \f$\int_D b_i\ \partial_y b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> y_grad_matrix;
    /// Stiffness matrix \f$\int_D \boldsymbol{\nabla} b_i\ \boldsymbol{\nabla} b_j\ dr\f$.
    Teuchos::SerialDenseMatrix<int,double> stiffness_matrix;
    /// Vector of the x component of the gradient matrices associated to edges.
    vector<Teuchos::SerialDenseMatrix<int,double> > x_grad_edge_matrix;
    /// Vector of the y component of the gradient matrices associated to edges.
    vector<Teuchos::SerialDenseMatrix<int,double> > y_grad_edge_matrix;
};

inline unsigned int FINITE_ELEMENT::Get_dof_per_cell() const
{
  return dof_per_cell;
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_downwind_matrix(
    unsigned int i) const
{
  return &downwind[i];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_upwind_matrix(
    unsigned int i) const
{
  return &upwind[i];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_edge_deln_matrix(
    unsigned int i,unsigned int j) const
{
  return &edge_deln_matrix[i][j];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_coupling_edge_deln_matrix(
    unsigned int i,unsigned int j) const
{
  return &coupling_edge_deln_matrix[i][j];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_mass_matrix() const
{
  return &mass_matrix;
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_x_grad_matrix() const
{
  return &x_grad_matrix;
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_y_grad_matrix() const
{
  return &y_grad_matrix;
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_x_grad_matrix(unsigned int i) const
{
  return &x_grad_edge_matrix[i];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_y_grad_matrix(unsigned int i) const
{
  return &y_grad_edge_matrix[i];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const FINITE_ELEMENT::Get_stiffness_matrix() const
{
  return &stiffness_matrix;
}

#endif
