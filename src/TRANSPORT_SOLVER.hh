#ifndef _TRANSPORT_SOLVER_HH_
#define _TRANSPORT_SOLVER_HH_

#include <ctime>
#include <fstream>
#include <string>
#include "AztecOO.h"
#include "Epetra_FEVector.h"
#include "Eper"
#include "Epetra_Map.h"
#include "Epetra_SerialComm.h"
#include "CELL.hh"
#include "DOF_HANDLER.hh"
#include "FINITE_ELEMENT.hh"
#include "GLC.hh"
#include "LS.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
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
    TRANSPORT_SOLVER(string* g_inputfile,string* p_inputfile,string* outputfile);

    ~TRANSPORT_SOLVER();

    /// Solve the transport equation.
    void Solve();

    /// Write the mesh and the solution in a file. Assume linear finite
    /// elements.
    void Write_in_file();

  private :
    /// Starting time of the calculation.
    time_t start_time;
    /// Ending time of the calculation.
    time_t end_time;
    /// Geometry input filename.
    string* geometry_inputfile;
    /// Parameters input filename.
    string* parameters_inputfile;
    /// Output filename.
    string* outputfile;
    /// Flux moments vector.
    Epetra_FeVector* flux_moments;
    /// Epetra communicator.
    Epetra_SerialComm comm;
    /// Epetra map associated to the #flux_moments.
    Epetra_Map* flux_moments_map;
    /// Degrees of freedom handler of the problem.
    DOF_HANDLER* dof_handler;
    /// Parameters of the problem.
    PARAMETERS parameters;
    /// Triangulation of the problems.
    TRIANGULATION triangulation;
    /// Vector of the different quadratures.
    vector<QUADRATURE*> quad;
};

#endif
