#include <cassert>
#include <string>
#include <vector>
#include "TRIANGULATION.hh"

using namespace std;

int main(int argc,char** argv)
{
  const unsigned int n_cells(3);
  const unsigned int n_mat(3);
  const unsigned int n_src(3);
  const CELL_TYPE cell_type(polygon);
  const double x_0(0.);
  const double x_1(0.5);
  const double x_2(1.);
  const double x_3(2.);
  const double y_0(0.);
  const double y_1(0.5);
  const double y_2(1.);
  const double y_3(2.);
  vector<unsigned int> n_vertices(3,0.);
  n_vertices[0] = 3;
  n_vertices[1] = 4;
  n_vertices[2] = 5;
  vector<unsigned int> mat_id(3,0);
  vector<unsigned int> src_id(3,0);
  mat_id[0] = 0;
  mat_id[1] = 1;
  mat_id[2] = 2;
  src_id[0] = 0;
  src_id[1] = 2;
  src_id[2] = 1;

  string filename("geometry_pol.inp");

  TRIANGULATION triangulation(&filename);
  triangulation.Read_geometry();
  triangulation.Build_edges();

  // Check the number of cells 
  assert(n_cells==triangulation.Get_n_cells());

  // Check the number of materials
  assert(n_mat==triangulation.Get_n_materials());

  // Check the number of sources
  assert(n_src==triangulation.Get_n_sources());

  // Check the number of vertices per cell
  for (unsigned int i=0; i<n_cells; ++i)
    assert(n_vertices[i]==triangulation.Get_n_vertices(i));

  // Check the coordinates of the vertices of the first cell
  assert(x_0==triangulation.Get_x_coord(0,0));
  assert(y_0==triangulation.Get_y_coord(0,0));
  assert(x_2==triangulation.Get_x_coord(0,1));
  assert(y_0==triangulation.Get_y_coord(0,1));
  assert(x_0==triangulation.Get_x_coord(0,2));
  assert(y_2==triangulation.Get_y_coord(0,2));
  
  // Check the coordinates of the vertices of the first cell
  assert(x_0==triangulation.Get_x_coord(1,0));
  assert(y_2==triangulation.Get_y_coord(1,0));
  assert(x_1==triangulation.Get_x_coord(1,1));
  assert(y_1==triangulation.Get_y_coord(1,1));
  assert(x_2==triangulation.Get_x_coord(1,2));
  assert(y_3==triangulation.Get_y_coord(1,2));
  assert(x_0==triangulation.Get_x_coord(1,3));
  assert(y_3==triangulation.Get_y_coord(1,3));

  // Check the coordinates of the vertices of the first cell
  assert(x_2==triangulation.Get_x_coord(2,0));
  assert(y_0==triangulation.Get_y_coord(2,0));
  assert(x_3==triangulation.Get_x_coord(2,1));
  assert(y_0==triangulation.Get_y_coord(2,1));
  assert(x_3==triangulation.Get_x_coord(2,2));
  assert(y_3==triangulation.Get_y_coord(2,2));
  assert(x_2==triangulation.Get_x_coord(2,3));
  assert(y_3==triangulation.Get_y_coord(2,3));
  assert(y_1==triangulation.Get_y_coord(2,4));
  assert(y_1==triangulation.Get_y_coord(2,4));

  // Check mat_id
  for (unsigned int i=0; i<n_cells; ++i)
    assert(mat_id[i]==triangulation.Get_mat_id(i));

  // Check src_id
  for (unsigned int i=0; i<n_cells; ++i)
    assert(src_id[i]==triangulation.Get_src_id(i));

  // Check the type of cells
  assert(cell_type==triangulation.Get_cell_type());

  return 0;
}
