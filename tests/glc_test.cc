#include <cassert>
#include <cmath>
#include "gsl_math.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "GLC.hh"

using namespace std;

int main(int argc,char** argv)
{
  unsigned int n_dir(12);
  unsigned int n_mom(15);
  const double four_pi(4.*M_PI);
  d_vector omega(3,0.);
  omega[0] = 0.868846143426105;
  omega[1] = 0.35988785622265201;
  omega[2] = 0.33998104358485631;

  GLC quad(4,4,false);
  quad.Build_quadrature(four_pi);

  // Check the number of direction
  assert(n_dir==quad.Get_n_dir());

  // Check the number of momemts
  assert(n_mom==quad.Get_n_mom());

  // Check the degree of expension
  assert(quad.Get_l(0)==0);
  assert(quad.Get_l(1)==1);
  assert(quad.Get_l(2)==1);
  assert(quad.Get_l(3)==2);
  assert(quad.Get_l(4)==2);
  assert(quad.Get_l(5)==2);
  assert(quad.Get_l(6)==3);
  assert(quad.Get_l(7)==3);
  assert(quad.Get_l(8)==3);
  assert(quad.Get_l(9)==3);
  assert(quad.Get_l(10)==4);
  assert(quad.Get_l(11)==4);
  assert(quad.Get_l(12)==4);
  assert(quad.Get_l(13)==4);
  assert(quad.Get_l(14)==4);

  // Check omega 
  Teuchos::SerialDenseVector<int,double> const* const omega_ptr(quad.Get_omega(0));
  assert(fabs(omega[0]-(*omega_ptr)(0))<1e-12);
  assert(fabs(omega[1]-(*omega_ptr)(1))<1e-12);
  assert(fabs(omega[2]-(*omega_ptr)(2))<1e-12);
}
