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

#ifndef _LS_HH_
#define _LS_HH_

#include "EXCEPTION.hh"
#include "QUADRATURE.hh"

using namespace std;

/**
 * This class build the Level Symmetric quadrature.
 */

class LS : public QUADRATURE
{
  public :
    LS(unsigned int sn,unsigned int L_max,bool galerkin);

  private :
    /// Compute omega in one octant.
    void Build_octant();

    /// Compute the different omega given a set of directions.
    void Compute_omega(d_vector const &direction);
};

#endif
