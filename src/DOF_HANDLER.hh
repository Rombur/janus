#ifndef _DOF_HANDLER_HH_
#define _DOF_HANDLER_HH_

#include <list>
#include <map>
#include <set>
#include <vector>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "BLD.hh"
#include "CELL.hh"
#include "EDGE.hh"
#include "FINITE_ELEMENT.hh"
#include "PARAMETERS.hh"
#include "PWLD.hh"
#include "QUADRATURE.hh"
#include "TRIANGULATION.hh"

using namespace std;

typedef list<unsigned int> ui_list;
typedef map<unsigned int,unsigned int> ui_map;
typedef set<unsigned int> ui_set;
typedef vector<unsigned int> ui_vector;

/**
 * DOF_HANDLER creates and links the cells. Build the different sweep
 * ordering.
 */
class DOF_HANDLER
{
  public :
    DOF_HANDLER(TRIANGULATION* triang,PARAMETERS &param);

    ~DOF_HANDLER();

    /// Compute the ordering of the cell for the sweep associated to each /// direction for each QUADRATURE object.
    void Compute_sweep_ordering(vector<QUADRATURE*> &quad);

    /// Return true if the direction is the most normal for the bottom
    /// boundary.
    bool Is_most_normal_bottom(unsigned int lvl,unsigned int direction) const;

    /// Return true if the direction is the most normal for the right
    /// boundary.
    bool Is_most_normal_right(unsigned int lvl,unsigned int direction) const;

    /// Return true if the direction is the most normal for the top
    /// boundary.
    bool Is_most_normal_top(unsigned int lvl,unsigned int direction) const;

    /// Return true if the direction is the most normal for the left
    /// boundary.
    bool Is_most_normal_left(unsigned int lvl,unsigned int direction) const;

    /// Return the number of significant angular fluxes per direction.
    unsigned int Get_n_sf_per_dir() const;

    /// Return the number of degrees of freedom.
    unsigned int Get_n_dof() const;

    /// Return the number of cells.
    unsigned int Get_n_cells() const;

    /// Return the index of first basis function of a "reflective" cell given 
    /// the cell id.
    unsigned int Get_saf_map_dof(unsigned int cell_id);

    /// Return the "reflective" index of first basis function of a "reflective" 
    /// cell given the cell id.
    unsigned int Get_saf_map_reflective_dof(unsigned int cell_id);

    /// Return the sweep ordering for a given quadrature and a given
    /// direction.
    ui_vector const* const Get_sweep_order(unsigned int q,unsigned int direction) 
      const;

    /// Return a pointer to a cell.
    CELL* Get_cell(unsigned int i);

    /// Return the begin iterator of the mesh vector.
    vector<CELL*>::iterator Get_mesh_begin();

    /// Return the end iterator of the mesh vector.
    vector<CELL*>::iterator Get_mesh_end();

    /// Get an iterator to the beginning of #edges.
    vector<EDGE>::iterator Get_edges_begin();

    /// Get an iterator to the end of #edges.
    vector<EDGE>::iterator Get_edges_end();

  private :
    /// Compute the most normal directions among the quadrature.
    void Compute_most_normal_direction(QUADRATURE* quad,unsigned int lvl);

    /// Compute the external normal associated to each edge of each cell.
    void Compute_external_normal();

    /// Number of significant angular fluxes per direction.
    unsigned int n_sf_per_dir;
    /// Number of degrees of freedoms.
    unsigned int n_dof;
    /// Significant angular fluxes map of the dof: the key is the cell id and the 
    /// value is the index of the first basis function.
    ui_map saf_map_dof;
    /// Significant angular fluxes map of the "reflective" dof: the key is the cell 
    /// id  and the value is the reflective index of the first basis function.
    ui_map saf_map_reflective_dof;
    /// Pointer to the triangulation of the geometry.
    TRIANGULATION *triangulation;
    /// Most normal directions on the bottom boundary for each quadrature.
    vector<ui_set> most_n_bottom;
    /// Most normal directions on the right boundary for each quadrature.
    vector<ui_set> most_n_right;
    /// Most normal directions on the top boundary for each quadrature.
    vector<ui_set> most_n_top;
    /// Most normal directions on the left boundary for each quadrature.
    vector<ui_set> most_n_left;
    /// Vector of cells representing the mesh.
    vector<CELL*> mesh;
    /// Sweep ordering associated to the different quadratures.
    vector<vector<ui_vector> > sweep_order;
};

inline bool DOF_HANDLER::Is_most_normal_bottom(unsigned int lvl,
    unsigned int direction) const
{
  return most_n_bottom[lvl].count(direction);
}

inline bool DOF_HANDLER::Is_most_normal_right(unsigned int lvl,
    unsigned int direction) const
{
  return most_n_right[lvl].count(direction);
}

inline bool DOF_HANDLER::Is_most_normal_top(unsigned int lvl,
    unsigned int direction) const
{
  return most_n_top[lvl].count(direction);
}

inline bool DOF_HANDLER::Is_most_normal_left(unsigned int lvl,
    unsigned int direction) const
{
  return most_n_left[lvl].count(direction);
}

inline unsigned int DOF_HANDLER::Get_n_sf_per_dir() const
{
  return n_sf_per_dir;
}

inline unsigned int DOF_HANDLER::Get_n_dof() const
{
  return n_dof;
}

inline unsigned int DOF_HANDLER::Get_n_cells() const
{
  return triangulation->Get_n_cells();
}

inline unsigned int DOF_HANDLER::Get_saf_map_dof(unsigned int cell_id) 
{
  return saf_map_dof[cell_id];
}

inline unsigned int DOF_HANDLER::Get_saf_map_reflective_dof(unsigned int cell_id) 
{
  return saf_map_reflective_dof[cell_id];
}

inline ui_vector const* const DOF_HANDLER::Get_sweep_order(unsigned int q,
    unsigned int direction) const
{
  return &sweep_order[q][direction];
}
    
inline CELL* DOF_HANDLER::Get_cell(unsigned int i) 
{
  return mesh[i];
}

inline vector<CELL*>::iterator DOF_HANDLER::Get_mesh_begin() 
{
  return mesh.begin();
}

inline vector<CELL*>::iterator DOF_HANDLER::Get_mesh_end() 
{
  return mesh.end();
}

inline vector<EDGE>::iterator DOF_HANDLER::Get_edges_begin()
{
  return triangulation->Get_edges_begin();
}

inline vector<EDGE>::iterator DOF_HANDLER::Get_edges_end()
{
  return triangulation->Get_edges_end();
}

#endif
