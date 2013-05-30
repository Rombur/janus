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

#include "GLC.hh"

GLC::GLC(unsigned int sn,unsigned int L_max,bool galerkin) :
  QUADRATURE(sn,L_max,galerkin)
{}

void GLC::Build_octant()
{
  d_vector azimuthal_nodes((sn*(sn+2))/8);
  d_vector azimuthal_weight((sn*(sn+2))/8);
  d_vector polar_weight(sn/2,0.);
  d_vector cos_theta(sn/2,0.);
  d_vector sin_theta(sn/2,0.);

  // This function determines the Gauss-Legendre abscissae and weights
  // necessary for an n-point fixed order integration scheme.
  gsl_integration_glfixed_table* table(gsl_integration_glfixed_table_alloc(sn));

  for (unsigned int i=sn/2; i<sn; ++i)
  {
    // This function applies the Gauss-Legendre integration rule contained in
    // the table and returns the result.
    gsl_integration_glfixed_point(-1.,1.,i,&cos_theta[i-sn/2],&polar_weight[i-sn/2],
       table);
    sin_theta[i-sn/2] = sqrt(1.-pow(cos_theta[i-sn/2],2));
  }

  // Build the Chebyshev quadrature
  Build_Chebyshev_quadrature(azimuthal_nodes,azimuthal_weight);

  unsigned int pos(0);
  unsigned int offset(0);
  for (unsigned int i=0; i<sn/2; ++i)
  {
    for (unsigned int j=0; j<sn/2-i; ++j)
    {
      omega[pos](0) = sin_theta[i]*cos(azimuthal_nodes[j+offset]);
      omega[pos](1) = sin_theta[i]*sin(azimuthal_nodes[j+offset]);
      omega[pos](2) = cos_theta[i];
      weight(pos) = polar_weight[i]*azimuthal_weight[j+offset];
      ++pos;
    }
    offset += sn/2-i;
  }

  // This function frees the memory associated with the table.
  gsl_integration_glfixed_table_free(table);
}

void GLC::Build_Chebyshev_quadrature(d_vector &nodes,d_vector &weight)
{
  unsigned int pos(0);
  for (unsigned int i=0; i<sn/2; ++i)
  {
    const double j_max(sn/2-i);
    for (double j=0; j<j_max; ++j)
    {
      nodes[pos] = M_PI_2*j/j_max+M_PI_4/j_max;
      weight[pos] = M_PI_2/j_max;
      ++pos;
    }
  }
}
