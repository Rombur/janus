#ifndef _TRIANGULATION_HH_
#define _TRIANGULATION_HH_

#include <fstream>
#include <string>
#include <vector>
#include "EDGE.hh"

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;

enum CELL_TYPE{rectangle,polygon};

class TRIANGULATION
{
  public :
    TRIANGULATION(string* geometry_inputfile);
    
    /// Read the geometry input file.
    void Read_geometry();

    /// Construct all the edges of the triangulation.
    void Build_edges();

  private :
    /// Get the EDDGE_TYPE of an edge, when its vertices are known.
    EDGE_TYPE Get_edge_type(d_vector const &vertex_0,d_vector const &vertex_1);

    /// Number of cells.
    unsigned int n_cells;
    /// Dimension of the problem. For now, only 2D problems are allowed.
    const unsigned int dim;
    /// Type of the cells: rectangle or polygon.
    CELL_TYPE cell_type;
    /// The domain is assumed to be a rectangle and bottom is y_min.
    double bottom_value;
    /// The domain is assumed to be a rectangle and right is x_max.
    double right_value;
    /// The domain is assumed to be a rectangle and top is y_max.
    double top_value;
    /// The domain is assumed to be a rectangle and right is x_min.
    double left_value;
    /// Pointer to the name of the input file containing the geometry of the
    /// problem.
    string* geometry_filename;
    /// Number of vertices of each cell.
    ui_vector n_vertices;
    /// Material id of each cell.
    ui_vector mat_id;
    /// Source id of each cell.
    ui_vector src_id;
    /// Vector containing all the edge.
    vector<EDGE> edges;
    /// grid contains the coordinates of each vertices for each cells. The
    /// data is stored by [cell,vertex,componant].
    vector<vector<d_vector> > grid;
};

#endif
