#include <cassert>
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

  edge_0.Set_reflective_boundary();
  edge_0.Set_lid(0,2);
  edge_1.Set_lid(0,2);
  edge_1.Set_lid(1,1);

  edge_0.Set_external_normal(0,normal_0);
  edge_1.Set_external_normal(0,normal_1);
  edge_1.Set_external_normal(1,normal_2);

  return 0;
}
