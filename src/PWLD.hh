#ifndef _PWLD_HH_
#define _PWLD_HH_

#include <cmath>
#include <vector>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "EDGE.hh"
#include "CELL.hh"
#include "FINITE_ELEMENT.hh"

using namespace std;

typedef vector<double> d_vector;

/**
 * This class build the PieceWise Linear Discontinuous finite elements.
 */

class PWLD : public FINITE_ELEMENT
{
  public :
    PWLD(d_vector const &cell_x,d_vector const &cell_y);

    /// Build the 1-dimensional finite elements (edge).
    void Build_fe_1d();

    /// Build the upwindind matrices. upwind_dof_per_cell is not used in BLD.
    void Build_upwind_matrices(CELL* cell,vector<CELL*> const &mesh);

    /// Build the 2-dimensional finite elements (cell).
    void Build_fe_2d();

  private :
    /// Abscissa of the center of the element.
    double x_c;
    /// Ordinate of the center of the element.
    double y_c;
    /// Store temporary the edge mass matrix on a sub-triangle.
    vector<vector<d_vector> > edge_mass_matrix;
    /// Store temporary the edge deln matrix on a sub-triangle.
    vector<vector<vector<d_vector> > > deln_matrix;
};

#endif
