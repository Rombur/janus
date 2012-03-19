#include <cassert>
#include <string>
#include "TRIANGULATION.hh"

using namespace std;

int main(int argc,char** argv)
{
  const unsigned int n_cells(4);
  const unsigned int n_vertices(4);
  const unsigned int n_mat(4);
  const unsigned int n_src(4);
  const double x_0(0.);
  const double x_1(1.);
  const double x_2(2.);
  const double y_0(0.);
  const double y_1(1.1);
  const double y_2(2.2);
  vector<unsigned int> mat_id(4,0);
  vector<unsigned int> src_id(4,0);
  mat_id[0] = 0;
  mat_id[1] = 0;
  mat_id[2] = 1;
  mat_id[3] = 2;
  src_id[0] = 2;
  src_id[1] = 1;
  src_id[2] = 0;
  src_id[3] = 0;
  // Need to use AT_DATA instead of hard coding the path
  string filename("/home/bruno/Documents/Transport/janus/tests/geometry_rec.txt");

  TRIANGULATION triangulation(&filename);
  triangulation.Read_geometry();
  triangulation.Build_edges();

  // Check the number of cells 
  assert(n_cells==triangulation.Get_n_cells());

  // Check the number of vertices per cell
  for (unsigned int i=0; i<n_cells; ++i)
    assert(n_vertices==triangulation.Get_n_vertices(i));

  // Check the number of materials
  assert(n_mat==triangulation.Get_n_materials());

  // Check the number of sources
  assert(n_src==triangulation.Get_n_sources());

  // Check the abscissae of the vertices
  assert(x_0==triangulation.Get_x_coord(0,0));
  assert(x_0==triangulation.Get_x_coord(0,3));
  assert(x_0==triangulation.Get_x_coord(2,0));
  assert(x_0==triangulation.Get_x_coord(2,3));

  assert(x_1==triangulation.Get_x_coord(0,1));
  assert(x_1==triangulation.Get_x_coord(0,2));
  assert(x_1==triangulation.Get_x_coord(1,0));
  assert(x_1==triangulation.Get_x_coord(1,3));
  assert(x_1==triangulation.Get_x_coord(2,1));
  assert(x_1==triangulation.Get_x_coord(2,2));
  assert(x_1==triangulation.Get_x_coord(3,3));
  assert(x_1==triangulation.Get_x_coord(3,0));

  assert(x_2==triangulation.Get_x_coord(1,1));
  assert(x_2==triangulation.Get_x_coord(1,2));
  assert(x_2==triangulation.Get_x_coord(3,1));
  assert(x_2==triangulation.Get_x_coord(3,2));

  // Check the ordinated of the vertices
  assert(y_0==triangulation.Get_y_coord(0,0));
  assert(y_0==triangulation.Get_y_coord(0,1));
  assert(y_0==triangulation.Get_y_coord(1,0));
  assert(y_0==triangulation.Get_y_coord(1,1));

  assert(y_1==triangulation.Get_y_coord(0,2));
  assert(y_1==triangulation.Get_y_coord(0,3));
  assert(y_1==triangulation.Get_y_coord(1,2));
  assert(y_1==triangulation.Get_y_coord(1,3));
  assert(y_1==triangulation.Get_y_coord(2,0));
  assert(y_1==triangulation.Get_y_coord(2,1));
  assert(y_1==triangulation.Get_y_coord(3,0));
  assert(y_1==triangulation.Get_y_coord(3,1));

  assert(y_2==triangulation.Get_y_coord(2,2));
  assert(y_2==triangulation.Get_y_coord(2,3));
  assert(y_2==triangulation.Get_y_coord(3,2));
  assert(y_2==triangulation.Get_y_coord(3,3));

  // Check mat_id
  for (unsigned int i=0; i<n_cells; ++i)
    assert(mat_id[i]==triangulation.Get_mat_id(i));

  // Check src_id
  for (unsigned int i=0; i<n_cells; ++i)
    assert(src_id[i]==triangulation.Get_src_id(i));

  return 0;
}
