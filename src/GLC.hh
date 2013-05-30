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
