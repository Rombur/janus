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

#ifndef _ERROR_ESTIMATOR_HH_
#define _ERROR_ESTIMATOR_HH_

#include <algorithm>
#include <cmath>
#include <list>
#include <map>
#include <vector>
#include "Epetra_MultiVector.h"
#include "DOF_HANDLER.hh"
#include "PARAMETERS.hh"

using std::fabs;
using std::list;
using std::map;
using std::max_element;
using std::pair;
using std::pow;
using std::vector;

typedef list<unsigned int> ui_list;
typedef vector<double> d_vector;
typedef vector<unsigned int> ui_vector;

namespace ERROR_ESTIMATOR
{
  /// Compute the error estimate and flag the cells.
  void Compute_refinement(const unsigned int n_groups,DOF_HANDLER *dof_handler,
      PARAMETERS const *const parameters,Epetra_MultiVector const* const group_flux,
      ui_set &cells_to_refine,ui_set &adjacent_cells,
      map<unsigned int,vector<vector<d_vector> > > &edge_to_refine);
}

#endif
