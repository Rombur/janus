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

#ifndef _TRIANGULATION_HH_
#define _TRIANGULATION_HH_

#include <algorithm>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "EDGE.hh"
#include "EXCEPTION.hh"

using namespace std;

typedef multimap<double,int> di_multimap;
typedef pair<double,int> di_pair;
typedef set<unsigned int> ui_set;
typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;
typedef vector<int> i_vector;
typedef vector<unsigned int> ui_vector;

/// Type of the cell: rectangle or arbitrary polygon.
enum CELL_TYPE{rectangle,polygon};

class TRIANGULATION
{
  public :
    TRIANGULATION(string* geometry_inputfile);
    
    /// Read the geometry input file.
    void Read_geometry();

    /// Construct all the edges of the triangulation.
    void Build_edges();

    /// Modify the grid by refining the cells whose cell IDs are in
    /// cells_to_refine. The cells whose IDs are in adjacent_cells have one
    /// or more of their edges refined. edge_to_refine contains the vertex
    /// which precedes (positive rotation) the new vertex. projection is the
    /// map of interpolation between the old and the new grid.
    /// @todo Should not be necessary to rebuild the edges.
    void Refine_mesh(ui_set const &cells_to_refine,ui_set const &adjacent_cells,
        map<unsigned int,vector<d_vector> > &edge_to_refine,
        vector<ui_vector> &projection){};
    void Refine_mesh(ui_set const &cells_to_refine,ui_set const &adjacent_cells,
        map<unsigned int,vector<vector<d_vector> > > &edge_to_refine,
        vector<ui_vector> &projection);

    /// Get the number of cells.
    unsigned int Get_n_cells() const;

    /// Get the number of cells along x when the rectangular cells are used.
    unsigned int Get_n_x_cells() const;

    /// Get the number of cells along y when the rectangular cells are used.
    unsigned int Get_n_y_cells() const;

    /// Get the number of vertices for a cell i.
    unsigned int Get_n_vertices(unsigned int i) const;

    /// Get the material id of of the cell i.
    unsigned int Get_mat_id(unsigned int i) const;

    /// Get the number of materials.
    unsigned int Get_n_materials() const;

    /// Get the source id of of the cell i.
    unsigned int Get_src_id(unsigned int i) const;

    /// Get the number of sources.
    unsigned int Get_n_sources() const;

    /// Get the type of the cells: rectangle or polygon.
    CELL_TYPE Get_cell_type() const;

    /// Get the abscissa of vertex j of the cell i.
    double Get_x_coord(unsigned int i,unsigned int j) const;

    /// Get the ordinate of vertex j of the cell i.
    double Get_y_coord(unsigned int i,unsigned int j) const;

    /// Get the global id of the edge j of the cell i.
    unsigned int Get_edge_gid(unsigned int i,unsigned int j) const;

    /// Get a pointer to an edge given the global id of the edge.
    EDGE* Get_edge(unsigned int i);

    /// Get an iterator to the beginning of #edges.
    vector<EDGE>::iterator Get_edges_begin();

    /// Get an iterator to the end of #edges.
    vector<EDGE>::iterator Get_edges_end();

  private :
    /// Get the #EDGE_TYPE of an edge, when its vertices are known.
    EDGE_TYPE Get_edge_type(d_vector const &vertex_0,d_vector const &vertex_1);

    void Project(unsigned int i,unsigned int dof,unsigned int &pos,
        unsigned int n_old_vertices,
        d_vector const &vertex,vector<ui_vector> &projection,
        vector<vector<d_vector> > const &coarsest_vertices);

    /// Number of cells.
    unsigned int n_cells;
    /// Number of cells along x when rectangular cells are used.
    unsigned int n_x_cells;
    /// Number of cells along y when rectangular cells are used.
    unsigned int n_y_cells;
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
    /// Contains the edge global id associated to each cell.
    vector<d_vector> cell_to_edge_gid;
    /// grid contains the coordinates of each vertices for each cells. The
    /// data is stored by [cell,vertex,component]. The first component is the
    /// abcissa, the second component is the ordinate, and the third component
    /// is coarse flag (0 if the vertex is coarse, 1 otherwise).
    vector<vector<d_vector> > grid;
};

inline unsigned int TRIANGULATION::Get_n_cells() const
{
  return n_cells;
}

inline unsigned int TRIANGULATION::Get_n_x_cells() const
{
  return n_x_cells;
}

inline unsigned int TRIANGULATION::Get_n_y_cells() const
{
  return n_y_cells;
}

inline unsigned int TRIANGULATION::Get_n_vertices(unsigned int i) const
{
  return n_vertices[i];
}

inline unsigned int TRIANGULATION::Get_mat_id(unsigned int i) const
{
  return mat_id[i];
}

inline unsigned int TRIANGULATION::Get_src_id(unsigned int i) const
{
  return src_id[i];
}

inline CELL_TYPE TRIANGULATION::Get_cell_type() const
{
  return cell_type;
}

inline double TRIANGULATION::Get_x_coord(unsigned int i,unsigned int j) const
{
  return grid[i][j][0];
}

inline double TRIANGULATION::Get_y_coord(unsigned int i,unsigned int j) const
{
  return grid[i][j][1];
}

inline unsigned int TRIANGULATION::Get_edge_gid(unsigned int i,unsigned int j) const
{
  return cell_to_edge_gid[i][j];
}

inline EDGE* TRIANGULATION::Get_edge(unsigned int i)
{
  return &edges[i];
}

inline vector<EDGE>::iterator TRIANGULATION::Get_edges_begin()
{
  return edges.begin();
}

inline vector<EDGE>::iterator TRIANGULATION::Get_edges_end()
{
  return edges.end();
}

#endif
