#ifndef _QUADRATURE_HH_
#define _QUADRATURE_HH_

#include <cmath>
#include <vector>
#include "gsl_math.h"
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
    void Build_quadrature();

  protected :
    /// Purely virtual function. Compute omega in one octant.
    virtual void Build_octant() = 0;

    /// Deploy the octant.
    void Deploy_octant();

    /// Compute the spherical harmonics and build the matrix M2D.
    void Compute_harmonics();

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
    /// Weights of the quadrature when a non-Galerkin quadrature is used.
    Teuchos::SerialDenseVector<int,double> weight;
    /// Moments to directions matrix.
    Teuchos::SerialDenseMatrix<int,double> M2D;
    /// Directions to moments matrix.
    Teuchos::SerialDenseMatrix<int,double> D2M;
    /// Vector of omega for each direction.
    vector<Teuchos::SerialDenseVector<int,double> > omega;
};

#endif
