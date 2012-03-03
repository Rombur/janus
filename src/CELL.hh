#ifndef _CELL_HH_
#define _CELL_HH_

#include <cmath>
#include <vector>
#include "gsl_math.h"
#include "EDGE.hh"

using namespace std;

typedef vector<unsigned int> ui_vector;
typedef vector<double> d_vector;

class FINITE_ELEMENT;
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
        unsigned int last_dof,double source,double sigma_t,d_vector sigma_s,
        vector<EDGE*> edges, FINITE_ELEMENT* fe);

    /// Return the cell id.
    unsigned int Get_id() const;

    /// Return the begin iterator of the cell_edge vector.
    vector<EDGE*>::iterator Get_cell_edges_begin();

    /// Return the end iterator of the cell_edge vector.
    vector<EDGE*>::iterator Get_cell_edges_end();

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
    bool Is_same_vertex(d_vector const &v1,d_vector const &v2);

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
    /// Intensity of the source in the cell.
    double source;
    /// Diffusion coefficient in the cell.
    double D;
    /// Total cross section in the cell.
    double sigma_t;
    /// Scattering cross section in the cell.
    d_vector sigma_s;
    /// Ortogonal length of the cell associated to each edge.
    d_vector orthogonal_length;
    /// Vector of the pointer to edges which compose the cell.
    vector<EDGE*> cell_edges;
    /// Pointer to the finite element associated to the cell.
    FINITE_ELEMENT* fe;
};

inline unsigned int CELL::Get_id() const
{
  return id;
}

inline vector<EDGE*>::iterator CELL::Get_cell_edges_begin() 
{
  return cell_edges.begin();
}

inline vector<EDGE*>::iterator CELL::Get_cell_edges_end()
{
  return cell_edges.begin();
}

inline FINITE_ELEMENT const* const CELL::Get_fe() const
{
  return fe;
}

#endif
