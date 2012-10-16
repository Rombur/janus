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

#ifndef _QUADRATURE_HH_
#define _QUADRATURE_HH_

#include <cmath>
#include <vector>
#include "gsl_sf_legendre.h"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

using namespace std;

typedef vector<double> d_vector;

class QUADRATURE
{
  public :
    QUADRATURE(unsigned int sn,unsigned int L_max,bool galerkin);

    /// Build the quadrature, i.e. M,D and omega (direction vector).
    void Build_quadrature(const double weight_sum);

    /// Return the number of directions of the quadrature.
    unsigned int Get_n_dir() const;

    /// Return the number of moments.
    unsigned int Get_n_mom() const;

    /// Return the degree l of expension coefficient of the scattering
    /// cross section given the number i of a moment.
    unsigned int Get_l(const unsigned int i) const;

    /// Return a pointer to omega for direction idir.
    Teuchos::SerialDenseVector<int,double> const* const Get_omega(unsigned int idir) 
      const;

    /// Return the first two components of omega for direction idir.
    Teuchos::SerialDenseVector<int,double> Get_omega_2d(unsigned int idir) const;

    /// Return the moment-to-discrete matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_M2D() const;

    /// Return the discrete-to-moment matrix.
    Teuchos::SerialDenseMatrix<int,double> const* const Get_D2M() const;

  protected :
    /// Purely virtual function. Compute omega in one octant.
    virtual void Build_octant() = 0;

    /// Deploy the octant.
    void Deploy_octant();

    /// Compute the spherical harmonics and build the matrix M2D.
    void Compute_harmonics(const double weight_sum);

    /// If flag is true, the quadrature is a Galerkin quadrature.
    const bool galerkin;
    /// Sn order of the quadrature.
    const unsigned int sn;
    /// L_max of the quadrature.
    const unsigned int L_max;
    /// Number of directions.
    unsigned int n_dir;
    /// Number of moments.
    unsigned int n_mom;
    /// Store the degree of the scattering cross section expansion associated
    /// to the number of the moment.
    d_vector moment_to_order;
    /// Weights of the quadrature when a non-Galerkin quadrature is used.
    Teuchos::SerialDenseVector<int,double> weight;
    /// Moments to directions matrix.
    Teuchos::SerialDenseMatrix<int,double> M2D;
    /// Directions to moments matrix.
    Teuchos::SerialDenseMatrix<int,double> D2M;
    /// Vector of omega for each direction.
    vector<Teuchos::SerialDenseVector<int,double> > omega;
};

inline unsigned int QUADRATURE::Get_n_dir() const
{
  return n_dir;
}

inline unsigned int QUADRATURE::Get_n_mom() const
{
  return n_mom;
}

inline unsigned int QUADRATURE::Get_l(const unsigned int i) const
{
  return moment_to_order[i];
}

inline Teuchos::SerialDenseVector<int,double> const* const QUADRATURE::Get_omega(unsigned int idir) const
{
  return &omega[idir];
}

inline Teuchos::SerialDenseMatrix<int,double> const* const QUADRATURE::Get_M2D() const
{
  return &M2D;
}

inline Teuchos::SerialDenseMatrix<int,double> const* const QUADRATURE::Get_D2M() const
{
  return &D2M;
}

#endif
