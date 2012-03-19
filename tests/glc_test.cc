#include <cassert>
#include <cmath>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "GLC.hh"

using namespace std;

int main(int argc,char** argv)
{
  unsigned int n_dir(12);
  d_vector omega(3,0.);
  omega[0] = 0.868846143426105;
  omega[1] = 0.35988785622265201;
  omega[2] = 0.33998104358485631;

  GLC quad(4,4,false);
  quad.Build_quadrature();

  // Check the number of direction
  assert(n_dir==quad.Get_n_dir());

  // Check omega 
  Teuchos::SerialDenseVector<int,double> const* const omega_ptr(quad.Get_omega(0));
  assert(fabs(omega[0]-(*omega_ptr)(0))<1e-12);
  assert(fabs(omega[1]-(*omega_ptr)(1))<1e-12);
  assert(fabs(omega[2]-(*omega_ptr)(2))<1e-12);
}
