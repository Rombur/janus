#ifndef _BLD_HH_
#define _BLD_HH_

#include <vector>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "CELL.hh"
#include "FINITE_ELEMENT.hh"

using namespace std;

typedef vector<double> d_vector;

/**
 * This class build the BiLinear Discontinuous finite elements. This class
 * should be used only with rectangular cells. If the mesh is built using the
 * polygon option, it is important that the first point of the cell is the
 * bottom left point.
 */

class BLD : public FINITE_ELEMENT
{
  public :
    BLD(d_vector const &cell_x,d_vector const &cell_y);

    /// Build the 1-dimensional finite elements (edge).
    void Build_fe_1d();

    /// Build the upwindind matrices. cell and mesh are not used in BLD.
    void Build_upwind_matrices(CELL* cell,vector<CELL*> const &mesh);

    /// Build the 2-dimensional finite elements (cell).
    void Build_fe_2d();

  private :
    /// Length of the horizontal sides.
    double delta_x;
    /// Length of the vertical sides.
    double delta_y;
};

#endif
