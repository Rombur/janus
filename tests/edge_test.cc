#include <cassert>
#include <cmath>
#include <vector>
#include "Teuchos_SerialDenseVector.hpp"
#include "EDGE.hh"

using namespace std;

typedef vector<double> d_vector;

int main(int argc,char** argv)
{
  const double v_0_x(0.303);
  const double v_0_y(0.503);
  const double v_1_x(-0.4);
  const double v_1_y(0.6);
  const double v_2_x(1.4);
  const double v_2_y(1.6);
  d_vector vertex_0(2,0.);
  d_vector vertex_1(2,0.);
  d_vector vertex_2(2,0.);
  Teuchos::SerialDenseVector<int,double> normal_0(2);
  Teuchos::SerialDenseVector<int,double> normal_1(2);
  Teuchos::SerialDenseVector<int,double> normal_2(2);

  vertex_0[0] = v_0_x;
  vertex_0[1] = v_0_y;
  vertex_1[0] = v_1_x;
  vertex_1[1] = v_1_y;
  vertex_2[0] = v_2_x;
  vertex_2[1] = v_2_y;

  EDGE edge_0(bottom_boundary,0,vertex_0,vertex_1);
  EDGE edge_1(interior,1,vertex_1,vertex_2);

  normal_0(0) = 3.4;
  normal_0(1) = 4.3;
  normal_1(0) = 5.4;
  normal_1(1) = 3.6;
  normal_2(0) = -5.4;
  normal_2(1) = -3.6;

  edge_0.Set_lid(0,2);
  edge_1.Set_lid(0,2);
  edge_1.Set_lid(1,1);

  edge_0.Set_cell_index(0,0);
  edge_0.Set_cell_index(1,1);
  edge_1.Set_cell_index(0,1);
  edge_1.Set_cell_index(1,2);

  edge_0.Set_exterior_normal(0,normal_0);
  edge_1.Set_exterior_normal(0,normal_1);
  edge_1.Set_exterior_normal(1,normal_2);

  // Check interior edge
  assert(edge_0.Is_interior()==false);
  assert(edge_1.Is_interior()==true);

  // Check the coordinates comparison
  assert(edge_0.Has_same_coord(vertex_2,vertex_1)==false);
  assert(edge_0.Has_same_coord(vertex_1,vertex_0)==true);

  // Check the type of edge
  assert(edge_0.Get_edge_type()==bottom_boundary);
  assert(edge_1.Get_edge_type()==interior);

  // Check the global id of the edge
  assert(edge_0.Get_gid()==0);
  assert(edge_1.Get_gid()==1);

  // Check the local id
  assert(edge_0.Get_lid(0)==2);
  assert(edge_1.Get_lid(0)==2);
  assert(edge_1.Get_lid(1)==1);

  // Check the cell index
  assert(edge_0.Get_cell_index(0)==0);
  assert(edge_0.Get_cell_index(1)==1);
  assert(edge_1.Get_cell_index(0)==1);
  assert(edge_1.Get_cell_index(1)==2);

  // Check the length of the edge
  assert(fabs(edge_0.Get_length()-0.70966)<0.0001);
  assert(fabs(edge_1.Get_length()-2.05912)<0.0001);

  // Check the coordinates of vertex_0
  assert(edge_0.Get_v0_x()==0.303);
  assert(edge_0.Get_v0_y()==0.503);
  // Check the coordinates of vertex_1
  assert(edge_0.Get_v1_x()==-0.4);
  assert(edge_0.Get_v1_y()==0.6);
  // Check the coordinates of vertex_0
  assert(edge_1.Get_v0_x()==-0.4);
  assert(edge_1.Get_v0_y()==0.6);
  // Check the coordinates of vertex_1
  assert(edge_1.Get_v1_x()==1.4);
  assert(edge_1.Get_v1_y()==1.6);

  // Check the exterior normal component
  assert(edge_0.Get_exterior_normal_component(0,0)==3.4);
  assert(edge_0.Get_exterior_normal_component(0,1)==4.3);
  assert(edge_1.Get_exterior_normal_component(0,0)==5.4);
  assert(edge_1.Get_exterior_normal_component(0,1)==3.6);
  assert(edge_1.Get_exterior_normal_component(1,0)==-5.4);
  assert(edge_1.Get_exterior_normal_component(1,1)==-3.6);

  return 0;
}
