#include <cassert>
#include <iostream>
#include <vector>
#include "BLD.hh"
#include "CELL.hh"
#include "EDGE.hh"
#include "FINITE_ELEMENT.hh"

using namespace std;

typedef vector<double> d_vector;

int main(int argc,char** argv)
{
  vector<d_vector> reordered_vertices;
  d_vector vertex_0(2,0.);
  d_vector vertex_1(2,0.);
  d_vector vertex_2(2,0.);
  d_vector vertex_3(2,0.);
  d_vector vertex_4(2,0.);
  d_vector vertex_5(2,0.);
  vertex_0[0] = 0.;
  vertex_0[1] = 0.;
  vertex_1[0] = 1.;
  vertex_1[1] = 0.;
  vertex_2[0] = 1.;
  vertex_2[1] = 1.;
  vertex_3[0] = 0.75;
  vertex_3[1] = 1.25;
  vertex_4[0] = 0.;
  vertex_4[1] = 1.;
  vertex_5[0] = 0.25;
  vertex_5[1] = 1.25;


  EDGE edge_0(bottom_boundary,0,vertex_0,vertex_1);
  EDGE edge_1(right_boundary,1,vertex_1,vertex_2);
  EDGE edge_2(interior,2,vertex_2,vertex_3);
  EDGE edge_3(top_boundary,3,vertex_3,vertex_5);
  EDGE edge_4(interior,4,vertex_5,vertex_4);
  EDGE edge_5(left_boundary,5,vertex_4,vertex_0);

  unsigned int cell_id(1);
  unsigned int n_vertices(6);
  unsigned int first_dof(4);
  unsigned int last_dof(11);
  double source(10.);
  d_vector sigma_t(3,0.);
  sigma_t[0] = 10.;
  sigma_t[1] = 5.;
  sigma_t[2] = 2.5;
  vector<d_vector> sigma_s(3,d_vector(3,0.));
  sigma_s[0][0] = 10.;
  sigma_s[0][1] = 7.5;
  sigma_s[0][2] = 5.5;
  sigma_s[1][0] = 4.;
  sigma_s[1][1] = 3.5;
  sigma_s[1][2] = 2.5;
  sigma_s[2][0] = 2.;
  sigma_s[2][1] = 1.5;
  sigma_s[2][2] = 0.5;
  d_vector cell_x(2,0.);
  d_vector cell_y(2,0.);
  cell_x[0] = 0.;
  cell_x[1] = 1.;
  cell_y[0] = 0.;
  cell_y[1] = 1.;
  FINITE_ELEMENT* fe = new BLD(cell_x,cell_y);
  vector<EDGE*> edges;
  edges.push_back(&edge_0);
  edges.push_back(&edge_1);
  edges.push_back(&edge_2);
  edges.push_back(&edge_3);
  edges.push_back(&edge_4);
  edges.push_back(&edge_5);
  
  CELL cell(cell_id,n_vertices,first_dof,last_dof,source,sigma_t,sigma_s,edges,fe);

  // Check the cell id
  assert(cell_id==cell.Get_id());

  // Check the number of edges
  assert(n_vertices==cell.Get_n_edges());

  // Check the first dof
  assert(first_dof==cell.Get_first_dof());

  // Check the last dof
  assert(last_dof==cell.Get_last_dof());

  // Check the source
  assert(source==cell.Get_source());

  // Check sigma_t
  assert(cell.Get_sigma_t(0)==10.);
  assert(cell.Get_sigma_t(1)==5.);
  assert(cell.Get_sigma_t(2)==2.5);

  // Check sigma_s
  assert(cell.Get_sigma_s(0,0)==10.);
  assert(cell.Get_sigma_s(0,1)==7.5);
  assert(cell.Get_sigma_s(0,2)==5.5);
  assert(cell.Get_sigma_s(1,0)==4.); 
  assert(cell.Get_sigma_s(1,1)==3.5);
  assert(cell.Get_sigma_s(1,2)==2.5);
  assert(cell.Get_sigma_s(2,0)==2.); 
  assert(cell.Get_sigma_s(2,1)==1.5);
  assert(cell.Get_sigma_s(2,2)==0.5);

  // Check the edge iterator
  vector<EDGE*>::iterator edge(cell.Get_cell_edges_begin());
  vector<EDGE*>::iterator edge_end(cell.Get_cell_edges_end());
  assert((*edge)->Get_gid()==0);
  ++edge;
  assert((*edge)->Get_gid()==1);
  ++edge;
  assert((*edge)->Get_gid()==2);
  ++edge;
  assert((*edge)->Get_gid()==3);
  ++edge;
  assert((*edge)->Get_gid()==4);
  ++edge;
  assert((*edge)->Get_gid()==5);
  ++edge;
  assert(edge==edge_end);

  // Check the reordering of the vertices
  reordered_vertices = cell.Reorder_vertices();
  assert(reordered_vertices[0][0]==0.);
  assert(reordered_vertices[0][1]==0.);
  assert(reordered_vertices[1][0]==1.);
  assert(reordered_vertices[1][1]==0.);
  assert(reordered_vertices[2][0]==1.);
  assert(reordered_vertices[2][1]==1.);
  assert(reordered_vertices[3][0]==0.75);
  assert(reordered_vertices[3][1]==1.25);
  assert(reordered_vertices[4][0]==0.25);
  assert(reordered_vertices[4][1]==1.25);
  assert(reordered_vertices[5][0]==0.);
  assert(reordered_vertices[5][1]==1.);

  return 0;
}
