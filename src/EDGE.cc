#include "EDGE.hh"

EDGE::EDGE(EDGE_TYPE e_type,unsigned int edge_gid,d_vector v0,d_vector v1) :
  edge_type(e_type),
  gid(edge_gid),
  vertex_0(v0),
  vertex_1(v1)
{
  length = sqrt(pow(vertex_0[0]-vertex_1[0],2)+pow(vertex_0[1]-vertex_1[1],2));
  cell_indices.resize(2,0);
  lid.resize(2,0);
  external_normal.resize(2,Teuchos::SerialDenseVector<int,double> (2));
}

bool EDGE::Is_interior() const
{
  bool inside(false);
  if (edge_type==interior)
    inside = true;

  return inside;
}

bool EDGE::Has_same_coord(d_vector &v0,d_vector &v1) const
{
  bool same(false);
  if ((v0[0]==vertex_0[0] && v0[1]==vertex_0[1] && v1[0]==vertex_1[0] &&
      v1[1]==vertex_1[1]) || (v1[0]==vertex_0[0] && v1[1]==vertex_0[1] &&
        v0[0]==vertex_1[0] && v0[1]==vertex_1[1]))
    same = true;

  return same;
}
