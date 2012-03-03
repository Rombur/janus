#ifndef _GLC_HH_
#define _GLC_HH_

#include <cmath>
#include <vector>
#include "gsl_integration.h"
#include "gsl_math.h"
#include "QUADRATURE.hh"

using namespace std;

typedef vector<double> d_vector;

/**
 * This class build the triangular Gauss-Legendre-Chebyshev quadrature.
 */ 

class GLC : public QUADRATURE
{
  public : 
    GLC(unsigned int sn,unsigned int L_max,bool galerkin);

  private :
    /// Compute omega in one octant.
    void Build_octant();

    /// Build the Chebyshev quadrature.
    void Build_Chebyshev_quadrature(d_vector &nodes,d_vector &weight);
};

#endif
