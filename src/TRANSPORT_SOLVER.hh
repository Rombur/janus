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

#ifndef _TRANSPORT_SOLVER_HH_
#define _TRANSPORT_SOLVER_HH_

#include <cassert>
#include <cmath>
#include <fstream>
#include <string>
#include "AztecOO.h"
#include "Epetra_BLAS.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_Time.hpp"
#include "CELL.hh"
#include "CROSS_SECTIONS.hh"
#include "DOF_HANDLER.hh"
#include "EXCEPTION.hh"
#include "FINITE_ELEMENT.hh"
#include "GLC.hh"
#include "LS.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
#include "TRANSPORT_OPERATOR.hh"
#include "TRIANGULATION.hh"

using namespace std;

/**
 * Solve the transport equation with the following options: angular multigrid
 * preconditioning, MIP preconditioning with cg preconditioned by jacobi or
 * symmetric Gauss-Seidel or ML, SI, BiCGSTAB, GMRES or GMRES with estimation
 * of the extrem eigenvalues.
 */
class TRANSPORT_SOLVER
{
  public :
    TRANSPORT_SOLVER(string* g_inputfile,string* p_inputfile,string* xs_inputfile,
        string* outputfile,Epetra_MpiComm* mpi_comm);

    ~TRANSPORT_SOLVER();

    /// Solve the transport equation.
    void Solve();

    /// Write the mesh and the solution in a file. Assume linear finite
    /// elements.
    void Write_in_file();

  private :
    /// Compute the convergence over a supergroup or over all the groups.
    double Compute_convergence(Epetra_MultiVector const &flux,
        Epetra_MultiVector const &old_flux, const unsigned int n) const;

    /// Size of the flux moments vector.
    unsigned int flux_moments_size;
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
    /// MultiVector of the groups of flux moments.
    Epetra_MultiVector* group_flux;
    /// Epetra map associated to the flux moments.
    Epetra_Map* flux_moments_map;
    /// Cross sections of the problem.
    CROSS_SECTIONS cross_sections;
    /// Parameters of the problem.
    PARAMETERS parameters;
    /// Triangulation of the problems.
    TRIANGULATION triangulation;
    /// Vector of the different quadratures.
    vector<QUADRATURE*> quad;
    /// Degrees of freedom handler of the problem.
    DOF_HANDLER* dof_handler;
};

#endif
