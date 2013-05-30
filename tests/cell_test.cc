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
  unsigned int n_groups(3);
  unsigned int n_levels(2);
  unsigned int L_max(4);
  d_vector source(n_groups,0.);
  d_vector sigma_e(n_groups,0.);
  vector<d_vector> sigma_t(n_groups,d_vector(n_levels,0.));
  vector<vector<vector<d_vector> > > sigma_s(n_groups, 
      vector<vector<d_vector> >(n_groups,vector<d_vector> (n_levels,
          d_vector(L_max+1,0.))));

  for (unsigned int i=0; i<n_groups; ++i)
    source[i] = double(i);

  for (unsigned int i=0; i<n_groups; ++i)
    sigma_e[i] = double(i);

  for (unsigned int i=0; i<n_groups; ++i)
    for (unsigned int j=0; j<n_levels; ++j)
      sigma_t[i][j] = double(i+j);

  for (unsigned int i=0; i<n_groups; ++i)
    for (unsigned int j=0; j<n_groups; ++j)
      for (unsigned int k=0; k<n_levels; ++k)
        for (unsigned int l=0; l<L_max+1; ++l)
          sigma_s[i][j][k][l] = double(i+j+k+l);

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
  
  CELL cell(cell_id,n_vertices,first_dof,last_dof,source,sigma_t,sigma_e,sigma_s,
      edges,fe);

  // Check the cell id
  assert(cell_id==cell.Get_id());

  // Check the number of edges
  assert(n_vertices==cell.Get_n_edges());

  // Check the first dof
  assert(first_dof==cell.Get_first_dof());

  // Check the last dof
  assert(last_dof==cell.Get_last_dof());

  // Check the source
  for (unsigned int i=0; i<n_groups; ++i)
    assert(cell.Get_source(i)==double(i));

  // Check sigma_e
  for (unsigned int i=0; i<n_groups; ++i)
    assert(cell.Get_sigma_e(i)==double(i));

  // Check sigma_t
  for (unsigned int i=0; i<n_groups; ++i)
    for (unsigned int j=0; j<n_levels; ++j)
      assert(cell.Get_sigma_t(i,j)==double(i+j));

  // Check sigma_s
  for (unsigned int i=0; i<n_groups; ++i)
    for (unsigned int j=0; j<n_groups; ++j)
      for (unsigned int k=0; k<n_levels; ++k)
        for (unsigned int l=0; l<L_max+1; ++l)
          assert(cell.Get_sigma_s(i,j,k,l)==double(i+j+k+l));

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
