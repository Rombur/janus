#ifndef _CELL_HH_
#define _CELL_HH_

#include <cmath>
#include <vector>
#include "gsl_math.h"
#include "EDGE.hh"
#include "FINITE_ELEMENT.hh"

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;

//class FINITE_ELEMENT;
/**
 *  Cell contains all the data which define a cell. If the number of edges is
 *  even, the orthogonal length is 2 X apothem (inradius) => assume
 *  that the polygon is regular (area = (apothem X perimeter)/2). If the
 *  number of edges is odd, the orthogonal length is the sum of the apothem
 *  and the circumradius => assume that the polygon is regular \f$area =
 *  \frac{r^2\times N \times \sin(2\pi /N)}{2}\f$
 */
class CELL
{
  public :
    CELL(unsigned int cell_id,unsigned int n_vertices,unsigned int first_dof,
        unsigned int last_dof,d_vector source,vector<d_vector> sigma_t,
        vector<vector<vector<d_vector> > > sigma_s,vector<EDGE*> edges,
        FINITE_ELEMENT* fe);

    ~CELL();

    /// Return the cell id.
    unsigned int Get_id() const;

    /// Return the number edges/vertices of the cell.
    unsigned int Get_n_edges() const;

    /// Return the index of the "first" basis function associated to the cell.
    unsigned int Get_first_dof() const;

    /// Return the index of the "last" basis function associated to the cell
    ///+1.
    unsigned int Get_last_dof() const;

    /// Return the intensity of the source in the cell for a given energy
    /// group.
    double Get_source(unsigned int g) const;

    /// Return the diffusion coefficient in the cell for a given energy group.
    double Get_diffusion_coefficient(unsigned int group) const;

    /// Return the \f$\Sigma_t\f$ in the cell for a given energy group and 
    /// angular level.
    double Get_sigma_t(unsigned int group,unsigned int lvl) const;

    /// Return the \f$\Sigma_s\f$ in the cell for a given energy group, angular 
    /// level, and moment.
    double Get_sigma_s(unsigned int group,unsigned int group_p,unsigned int lvl,
        unsigned int mom) const;
          
    /// Return the orthogonal length of the cell associated to the ith edge of
    /// the cell.
    double Get_orthogonal_length(unsigned int i);
    
    /// Return the begin iterator of the cell_edge vector.
    vector<EDGE*>::iterator Get_cell_edges_begin();

    /// Return the end iterator of the cell_edge vector.
    vector<EDGE*>::iterator Get_cell_edges_end();

    /// Return a pointer to fe.
    FINITE_ELEMENT* Get_mutable_fe();

    /// Return a const pointer to the constant fe (pointer and object are
    /// constant);
    FINITE_ELEMENT const* const Get_fe() const;

    /// Reorder the vertices to compute the area of the cell and for the
    /// output.
    vector<d_vector> Reorder_vertices();

  private :
    /// Compute the area and the perimeter of the cell.
    void Compute_area_and_perimeter();

    /// Return true if the two vertices are the same.
    bool Is_same_vertex(d_vector const* v1,d_vector const* v2) const;

    /// Cell id.
    unsigned int id;
    /// Number of vertices of the cell.
    unsigned int n_vertices;
    /// Index of the "first" basis function associated to the cell.
    unsigned int first_dof;
    /// Index of the "last" basis function associated to the cell + 1.
    unsigned int last_dof;
    /// Perimeter of the cell.
    double perimeter;
    /// Area of the cell.
    double area;
    /// Intensity of the source in the cell. source is a vector because of
    /// energy groups.
    d_vector source;
    /// Diffusion coefficient in the cell. D is a vector because of energy
    /// groups.
    d_vector D;
    /// Ortogonal length of the cell associated to each edge.
    d_vector orthogonal_length;
    /// Vector of the pointer to edges which compose the cell.
    vector<EDGE*> cell_edges;
    /// Total cross section in the cell. sigma_t is a vector of vector because of 
    /// energy groups and angular multigrid.
    vector<d_vector> sigma_t;
    /// Scattering cross section in the cell. sigma_s is a vector of vector of
    /// vector because of energy groups, the angular multigrid, and the
    /// Legendre expansion of the scattering cross sections.
    vector<vector<vector<d_vector> > > sigma_s;
    /// Finite elements associated to the cell.
    FINITE_ELEMENT* fe;
};

inline unsigned int CELL::Get_id() const
{
  return id;
}

inline unsigned int CELL::Get_n_edges() const
{
  return n_vertices;
}

inline unsigned int CELL::Get_first_dof() const
{
  return first_dof;
}

inline unsigned int CELL::Get_last_dof() const
{
  return last_dof;
}

inline double CELL::Get_source(unsigned int group) const
{
  return source[group];
}

inline double CELL::Get_diffusion_coefficient(unsigned int group) const
{
  return D[group];
}

inline double CELL::Get_sigma_t(unsigned int group,unsigned int lvl) const
{
  return sigma_t[group][lvl];
}

inline double CELL::Get_sigma_s(unsigned int group,unsigned int group_p,
    unsigned int lvl,unsigned int mom) const
{
  return sigma_s[group][group_p][lvl][mom];
}
inline double CELL::Get_orthogonal_length(unsigned int i)
{
  return orthogonal_length[i];
}

inline vector<EDGE*>::iterator CELL::Get_cell_edges_begin() 
{
  return cell_edges.begin();
}

inline vector<EDGE*>::iterator CELL::Get_cell_edges_end()
{
  return cell_edges.end();
}

inline FINITE_ELEMENT* CELL::Get_mutable_fe()
{
  return fe;
}

inline FINITE_ELEMENT const* const CELL::Get_fe() const
{
  return fe;
}

#endif
