#include <fstream>
#include <string>
#include <vector>

#ifndef _TRIANGULATION_HH_
#define _TRIANGULATION_HH_

using namespace std;

typedef vector<double> d_vector;

enum CELL_TYPE{rectangle,polygon};

class TRIANGULATION
{
  public :
    TRIANGULATION(string &geometry_inputfile);
    /// Read the geometry input file.
    void Read_geometry();

  private :
    /// Number of cells.
    unsigned int n_cells;
    /// Type of the cells: rectangle or polygon.
    CELL_TYPE cell_type;
    /// Number of vertices of each cell.
    ui_vector n_vertices;
};

#endif
