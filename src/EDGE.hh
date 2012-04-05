#ifndef _EDGE_HH_
#define _EDGE_HH_

#include <algorithm>
#include <cmath>
#include <vector>
#include "Teuchos_SerialDenseVector.hpp"

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

    /// Set the flag to flag to true if the edge is on the reflective
    /// boundary.
    void Set_reflective_boundary();

    /// Return true if the edge is on the reflective boundary.
    bool Is_reflective() const;

    /// Return true if the coordinates are the same that the one of the edge.
    /// The order does not matter.
    bool Has_same_coord(d_vector &v0,d_vector &v1) const;

    /// Return the type of the edge type (interior, bottom_boundary,
    /// right_bloundary, top_boundary or left_boundary).
    EDGE_TYPE Get_edge_type() const;

    /// Get the global id of the edge.
    unsigned int Get_gid() const;

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

    /// Get the abscissa of #vertex_0.
    double Get_v0_x() const;

    /// Get the ordinate of #vertex_0.
    double Get_v0_y() const;

    /// Get the abscissa of #vertex_1.
    double Get_v1_x() const;

    /// Get the ordinate of #vertex_1.
    double Get_v1_y() const;

    /// Get a pointer to #vertex_0.
    d_vector const* const Get_v0() const;

    /// Get a pointer to #vertex_1.
    d_vector const* const Get_v1() const;

    /// Set the normal in the ith component of #external_normal.
    void Set_external_normal(unsigned int i,
        Teuchos::SerialDenseVector<int,double> const &normal);

    /// Get the normal in the ith component of #external_normal.
    Teuchos::SerialDenseVector<int,double> const* const Get_external_normal(unsigned int i) const;

    /// Get the jth component of the normal in the ith component of #external_normal.
    double Get_external_normal_component(unsigned int i,unsigned int j) const;

  private :
    /// Flag to know if the edge is inside the medium, on the bottom boundary, on the
    /// the right boundary, on the top boundary or on the left boundary.
    EDGE_TYPE type;
    /// Flag for the reflective boundary.
    bool reflective_boundary;
    /// Global ID of the edge.
    unsigned int gid;
    /// Length of the edge.
    double length;
    /// Local id of the edge.
    ui_vector lid;
    /// Ids of the two cell associated to the edge.
    ui_vector cell_indices;
    /// Coordinates of the first vertex of the edge.
    d_vector vertex_0;
    /// Coordinates of the second vertex of the edge.
    d_vector vertex_1;
    /// Contains the normals associated to the edge for the two cells
    /// associated to the edge.
    vector<Teuchos::SerialDenseVector<int,double> > external_normal;
};

inline void EDGE::Set_reflective_boundary()
{
  reflective_boundary = true;
}

inline bool EDGE::Is_reflective() const
{
  return reflective_boundary;
}

inline EDGE_TYPE EDGE::Get_edge_type() const
{
  return type;
}

inline unsigned int EDGE::Get_gid() const
{
  return gid;
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

inline double EDGE::Get_v0_x() const
{
  return vertex_0[0];
}

inline double EDGE::Get_v0_y() const
{
  return vertex_0[1];
}

inline double EDGE::Get_v1_x() const
{
  return vertex_1[0];
}

inline double EDGE::Get_v1_y() const
{
  return vertex_1[1];
}

inline d_vector const* const EDGE::Get_v0() const
{
  return &vertex_0;
}

inline d_vector const* const EDGE::Get_v1() const
{
  return &vertex_1;
}

inline void EDGE::Set_external_normal(unsigned int i,
    Teuchos::SerialDenseVector<int,double> const &norm)
{
  external_normal[i] = norm;
}

inline Teuchos::SerialDenseVector<int,double> const* const EDGE::Get_external_normal(unsigned int i) const
{
  return &external_normal[i];
}

inline double EDGE::Get_external_normal_component(unsigned int i,unsigned int j) const
{
  return external_normal[i](j);
}

#endif
