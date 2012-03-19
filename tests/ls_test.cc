#include <cassert>
#include <cmath>
#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "LS.hh"

using namespace std;

int main(int argc,char** argv)
{
  unsigned int n_dir(12);
  d_vector omega(3,0.);
  omega[0] = 0.868890300722201205229788;
  omega[1] = 0.350021174581540677777041;
  omega[2] = 0.350021174581540677777041;
  
  LS quad(4,4,true);
  quad.Build_quadrature();

  // Check the number of direction
  assert(n_dir==quad.Get_n_dir());

  // Check omega and omega_2d
  Teuchos::SerialDenseVector<int,double> const* const omega_ptr(quad.Get_omega(0));
  Teuchos::SerialDenseVector<int,double> omega_2d(quad.Get_omega_2d(0));
  assert(omega[0]==(*omega_ptr)(0));
  assert(omega[1]==(*omega_ptr)(1));
  assert(omega[2]==(*omega_ptr)(2));
  assert(omega[0]==omega_2d(0));
  assert(omega[1]==omega_2d(1));

  // Check Galerkin
  Teuchos::BLAS<int,double> blas;
  Teuchos::SerialDenseMatrix<int,double> result(n_dir,n_dir);
  Teuchos::SerialDenseMatrix<int,double> const* const M2D(quad.Get_M2D());
  Teuchos::SerialDenseMatrix<int,double> const* const D2M(quad.Get_D2M());
  blas.GEMM(Teuchos::NO_TRANS,Teuchos::NO_TRANS,n_dir,n_dir,n_dir,1.,M2D->values(),
      M2D->stride(),D2M->values(),D2M->stride(),0.,result.values(),result.stride());
  for (unsigned int i=0; i<n_dir; ++i)
    for (unsigned int j=0; j<n_dir; ++j)
    {
      if (i==j) 
        assert(fabs(result(i,j)-1.)<1e-12);
      else
        assert(fabs(result(i,j))<1e-12);
    }

  return 0;
}
