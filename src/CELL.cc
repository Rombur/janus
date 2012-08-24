#include "CELL.hh"

CELL::CELL(unsigned int cell_id,unsigned int n_vertices,unsigned int first_dof,
        unsigned int last_dof,double source,d_vector sigma_t,vector<d_vector> sigma_s,
        vector<EDGE*> edges,FINITE_ELEMENT* fe) :
  id(cell_id),
  n_vertices(n_vertices),
  first_dof(first_dof),
  last_dof(last_dof),
  source(source),
  sigma_t(sigma_t),
  sigma_s(sigma_s),
  cell_edges(edges),
  fe(fe)
{
  unsigned int n_groups(sigma_t.size());
  unsigned int n_level(sigma_t[0].size());
  D.resize(n_groups);
  for (unsigned int g=0; g<n_groups; ++g)
  {
    if (sigma_s[n_level-1].size()>1)
      D[g] = 1./(3.*(sigma_t[g][n_level-1]-sigma_s[g][g][n_level-1][1]));
    else
      D[g] = 1./(3.*sigma_t[g][n_level-1]);
  }

  // Compute the area and the perimeter of the cell.
  Compute_area_and_perimeter();

  // Compute the length of the cell in the direction orthogonal to the edge.
  orthogonal_length.resize(n_vertices);
  if (n_vertices==3)
    for (unsigned int i=0; i<n_vertices; ++i)
      orthogonal_length[i] = 2.*area/cell_edges[i]->Get_length();
  else
  {
    if (n_vertices==4)
      for (unsigned int i=0; i<n_vertices; ++i)
        orthogonal_length[i] = area/cell_edges[i]->Get_length();
    else
    {
      if (n_vertices%2==0)
        for (unsigned int i=0; i<n_vertices; ++i)
          orthogonal_length[i] = 4.*area/perimeter;
      else
        for(unsigned int i=0; i<n_vertices; ++i)
          orthogonal_length[i] = 2.*area/perimeter+sqrt(2.*area/(n_vertices*
                sin(2.*M_PI/n_vertices)));
    }
  }
}

CELL::~CELL()
{
  // Delete the memory allocated in the DOF_HANDLER
  delete fe;
  fe = NULL;
}

vector<d_vector> CELL::Reorder_vertices()
{
  unsigned int index(0);
  vector<d_vector const*> vertices_done;
  vector<d_vector> points(n_vertices,(d_vector(2,0.)));
  for (unsigned int i=0; i<n_vertices; ++i)
  {
    bool v0_done(false);
    bool v1_done(false);
    bool v0_switch(false);
    bool v1_switch(false);
    const unsigned int vertices_done_size(vertices_done.size());
    for (unsigned int j=0; j<vertices_done_size; ++j)
    {
      if (Is_same_vertex(cell_edges[i]->Get_v0(),vertices_done[j]))
      {
        v0_done = true;
        if (j!=vertices_done_size-1)
          v0_switch = true;
      }
      if (Is_same_vertex(cell_edges[i]->Get_v1(),vertices_done[j]))
      {
        v1_done = true;
        if (j!=vertices_done_size-1)
          v1_switch = true;
      }
    }
    if (((v0_done==true) && (v0_switch==true) && (v1_done==false)) ||
        ((v1_done==true) && (v1_switch==true) && (v0_done==false)))
    {
      double tmp_0(points[index-1][0]);
      double tmp_1(points[index-1][1]);
      points[index-1][0] = points[index-2][0];
      points[index-1][1] = points[index-2][1];
      points[index-2][0] = tmp_0;
      points[index-2][1] = tmp_1;
    }
    if (v0_done==false)
    {
      points[index] = *(cell_edges[i]->Get_v0());
      vertices_done.push_back(cell_edges[i]->Get_v0());
      ++index;
    }
    if (v1_done==false)
    {
      points[index] = *(cell_edges[i]->Get_v1());
      vertices_done.push_back(cell_edges[i]->Get_v1());
      ++index;
    }
  }

  return points;
}

void CELL::Compute_area_and_perimeter()
{
  area = 0.;
  vector<d_vector> reordered_vertices(Reorder_vertices());
  for (int i=n_vertices-1; i>-1; --i)
  {
    int j(i+1);
    if (i==int(n_vertices-1))
      j=0;
    area += (reordered_vertices[j][0]+reordered_vertices[i][0])*
      (reordered_vertices[j][1]-reordered_vertices[i][1]);
  }
  area /= 2.;

  perimeter = 0.;
  for (unsigned int i=0; i<n_vertices; ++i)
    perimeter += cell_edges[i]->Get_length();
}

bool CELL::Is_same_vertex(d_vector const* v1,d_vector const* v2) const
{
  bool same(false);
  if (((*v1)[0]==(*v2)[0]) && ((*v1)[1]==(*v2)[1]))
    same = true;

  return same;
}
