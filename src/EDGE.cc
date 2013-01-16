/*
Copyright (c) 2012, Bruno Turcksin.

This file is part of Janus.

Janus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
he Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Janus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Janus.  If not, see <http://www.gnu.org/licenses/>.
*/

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
  exterior_normal.resize(2,Teuchos::SerialDenseVector<int,double> (2));
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
