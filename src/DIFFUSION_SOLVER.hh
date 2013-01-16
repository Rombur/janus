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

#ifndef _DIFFUSION_SOLVER_HH_
#define _DIFFUSION_SOLVER_HH_

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_Time.hpp"
#include "CROSS_SECTIONS.hh"
#include "DOF_HANDLER.hh"
#include "ERROR_ESTIMATOR.hh"
#include "MIP.hh"
#include "PARAMETERS.hh"
#include "TRIANGULATION.hh"

using namespace std;

typedef set<unsigned int> ui_set;
typedef vector<double> d_vector;
typedef vector<unsigned int> ui_vector;

/**
 * Solve the multigroup diffusion equation using MIP with adaptive mesh
 * refinement.
 */
class DIFFUSION_SOLVER
{
  public :
    DIFFUSION_SOLVER(string* g_inputfile,string* p_inputfile,string* xs_input_file,
        string* outputfile,Epetra_MpiComm* mpi_comm);

    ~DIFFUSION_SOLVER();

    /// Solve the diffusion equation.
    void Solve();

    /// Write the mesh and the solution in a file. Assume linear finite
    /// elements.
    void Write_in_file(string* filename=NULL);

  private :
    /// Refine the mesh and project the solution on the new grid.
    void Refine_mesh(unsigned int r);

    /// Project the solution on the new grid.
    void Project_solution(vector<ui_vector> const &projection);

    /// Compute the convergence over all the groups.
    double Compute_convergence(Epetra_MultiVector const &flux,
        Epetra_MultiVector const &old_flux, const unsigned int n) const;

    /// Size of the scalar flux vector.
    unsigned int scalar_flux_size;
    /// Timer for the initialization.
    Teuchos::Time* init_timer;
    /// Timer for the calculation.
    Teuchos::Time* calc_timer;
    /// Geometry input filename.
    string* geometry_inputfile;
    /// Parameters input filename.
    string* parameters_inputfile;
    /// Output filename.
    string* outputfile;
    /// Epetra communicator.
    Epetra_MpiComm* comm;
    /// MultiVector of the groups of scalar fluxes.
    Epetra_MultiVector* group_flux;
    /// Epetra map associated to the scalar fluxes.
    Epetra_Map* scalar_flux_map;
    /// Cross sections of the problem.
    CROSS_SECTIONS cross_sections;
    /// Parameters of the problem.
    PARAMETERS parameters;
    /// Triangulation of the problems.
    TRIANGULATION triangulation;
    /// Degrees of freedom handler of the problem.
    DOF_HANDLER* dof_handler;
};

#endif
