#ifndef _EDGE_HH_
#define _EDGE_HH_

#include <cmath>
#include <vector>

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;

/// Enum on the different types of edges.
enum EDGE_TYPE{interior,bottom_boundary,right_boundary,top_boundary,left_boundary};

/**
 * Edge contains data to represent an edge of the mesh. Each edge is own by two cells.
 */
class EDGE
{
  public :
    EDGE(EDGE_TYPE edge_type,unsigned int edge_gid,d_vector v0,d_vector v1);
    
    /// Return true if the edge is interior of the domain.
    bool Is_interior() const;

    /// Return true if the edge is on the reflective boundary.
    bool Is_reflective() const;

    /// Return true if the coordinates are the same that the one of the edge.
    /// The order does not matter.
    bool Has_same_coord(d_vector &v0,d_vector &v1) const;

    /// Set the local id of the edge in the ith component of lid.
    void Set_lid(unsigned int i,unsigned int cell_id);

    /// Get the local id of the edge in the ith component of lid.
    unsigned int Get_lid(unsigned int i) const;

    /// Set the cell id in the ith component of cell_indices.
    void Set_cell_index(unsigned int i,unsigned int cell_id);

    /// Get the cell id in the ith component of cell indices.
    unsigned int Get_cell_index(unsigned int i) const;

    /// Get the length of the edge.
    double Get_length() const;

  private :
    /// Flag to know if the edge is inside the medium, on the bottom boundary, on the
    /// the right boundary, on the top boundary or on the left boundary.
    const EDGE_TYPE type;
    /// Flag for the reflective boundary.
    bool reflective_boundary;
    /// Global ID of the edge.
    const unsigned int gid;
    /// Length of the edge.
    double length;
    /// Local id of the edge.
    ui_vector lid;
    /// Ids of the two cell associated to the edge.
    ui_vector cell_indices;
    /// Coordinates of the first vertex of the edge.
    const d_vector vertex_0;
    /// Coordinates of the second vertex of the edge.
    const d_vector vertex_1;
    /// Contains the normals associated to the edge for the two cells
    /// associated to the edge.
    vector<d_vector> external_normal;
};

// inline functions must be placed in .hh

inline bool EDGE::Is_reflective() const
{
  return reflective_boundary;
}

inline void EDGE::Set_lid(unsigned int i,unsigned int local_id)
{
  lid[i] = local_id;
}

inline unsigned int EDGE::Get_lid(unsigned int i) const
{
  return lid[i];
}

inline void EDGE::Set_cell_index(unsigned int i,unsigned int cell_id)
{
  cell_indices[i] = cell_id;
}

inline unsigned int EDGE::Get_cell_index(unsigned int i) const
{
  return cell_indices[i];
}

inline double EDGE::Get_length() const
{
  return length;
}

#endif
