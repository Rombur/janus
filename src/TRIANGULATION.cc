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

#include "TRIANGULATION.hh"

TRIANGULATION::TRIANGULATION(string* geometry_inputfile) : 
  n_x_cells(0),
  n_y_cells(0),
  dim(2),
  bottom_value(1e100),
  right_value(-1e100),
  top_value(-1e100),
  left_value(1e100),
  geometry_filename(geometry_inputfile)
{}

void TRIANGULATION::Read_geometry()
{
  // Open the file to read it
  ifstream geometry_file(geometry_filename->c_str(),ios::in);

  // Check that the file was open properly
  Check(geometry_file.good(),string ("Unable to open the file " + 
        *geometry_filename + " containing the geometry."));

  // Read the cell type: quadrilateral or polygon
  string cell_type_str;
  geometry_file>>cell_type_str;

  if (cell_type_str.compare("rectangle")==0 || 
      cell_type_str.compare("RECTANGLE")==0)
    cell_type = rectangle;
  else
    if (cell_type_str.compare("polygon")==0 ||
        cell_type_str.compare("POLYGON")==0)
      cell_type = polygon;

  if (cell_type==rectangle)
  {
    // Read the number of cells in x and y
    geometry_file>>n_x_cells>>n_y_cells;
    n_cells = n_x_cells*n_y_cells;
    mat_id.resize(n_cells);
    src_id.resize(n_cells);
    // Fill n_vertices with 4 since all the cells are rectangular
    n_vertices.resize(n_cells,4);
    // Read the abscissae
    const unsigned int n_x_vertices(n_x_cells+1);
    d_vector x(n_x_vertices,0.);
    for (unsigned int i=0; i<n_x_vertices; ++i)
    {
      geometry_file>>x[i];
      if (x[i]<left_value)
        left_value = x[i];
      if (x[i]>right_value)
        right_value = x[i];
    }
    // Read the ordinates
    const unsigned int n_y_vertices(n_y_cells+1);
    d_vector y(n_y_vertices,0.);
    for (unsigned int i=0; i<n_y_vertices; ++i)
    {
      geometry_file>>y[i];
      if (y[i]<bottom_value)
        bottom_value = y[i];
      if (y[i]>top_value)
        top_value = y[i]; 
    }
    // Read the material ids
    for (unsigned int i=0; i<n_cells; ++i)
      geometry_file>>mat_id[i];
    // Read the source ids
    for (unsigned int i=0; i<n_cells; ++i)
      geometry_file>>src_id[i];
    // Fill in grid
    grid.resize(n_cells);
    unsigned int x_pos(0),y_pos(0);
    for (unsigned int i=0; i<n_cells; ++i) 
    {
      vector<d_vector> rectangle(4,d_vector(dim+1,0.));
      if (i%n_x_cells==0 && i!=0)
      {
        x_pos = 0;
        y_pos += 1;
      }
      // Bottom left vertex
      rectangle[0][0] = x[x_pos];
      rectangle[0][1] = y[y_pos];
      // Bottom right vertex
      rectangle[1][0] = x[x_pos+1];
      rectangle[1][1] = y[y_pos];
      // Top right vertex
      rectangle[2][0] = x[x_pos+1];
      rectangle[2][1] = y[y_pos+1];
      // Top left vertex
      rectangle[3][0] = x[x_pos];
      rectangle[3][1] = y[y_pos+1];
      // Add the cell to the grid
      grid[i] = rectangle;
      ++x_pos;
    }
  }
  else 
  {
    // Read the number of cells
    geometry_file>>n_cells;
    n_vertices.resize(n_cells,0);
    grid.resize(n_cells);
    mat_id.resize(n_cells);
    src_id.resize(n_cells);
    for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      unsigned int n_vert_polygon(0);
      geometry_file>>n_vert_polygon;
      n_vertices[cell] = n_vert_polygon;
      vector<d_vector> polygon(n_vert_polygon);
      for (unsigned int vert=0; vert<n_vert_polygon; ++vert)
      {
        d_vector coord(dim+1,0.);
        for (unsigned int i=0; i<dim; ++i)
          geometry_file>>coord[i];
        if (coord[0]<left_value)
          left_value = coord[0];
        if (coord[0]>right_value)
          right_value = coord[0];
        if (coord[1]<bottom_value)
          bottom_value = coord[1];
        if (coord[1]>top_value)
          top_value = coord[1];
        polygon[vert] = coord;
      }
      grid[cell] = polygon;
      // Read the material id
      geometry_file>>mat_id[cell];
      // Read the source id
      geometry_file>>src_id[cell];
    }
  }
  geometry_file.close();
}

void TRIANGULATION::Build_edges()
{
  unsigned int edge_gid(0);
  cell_to_edge_gid.resize(n_cells);
  di_multimap vx_to_edge_pos;
  di_multimap vy_to_edge_pos;
  // Loop over the cells.
  for (unsigned int i=0; i<n_cells; ++i)
  {
    // Loop over the vertices of the cell.
    for (unsigned int vertex=0; vertex<n_vertices[i]; ++vertex)
    {
      d_vector vertex_0(dim+1,0.);
      d_vector vertex_1(dim+1,0.);
      if (vertex!=n_vertices[i]-1)
      {
        vertex_0 = grid[i][vertex];
        vertex_1 = grid[i][vertex+1];
      }
      else
      {
        vertex_0 = grid[i][vertex];
        vertex_1 = grid[i][0];
      }
      bool edge_exist(false);
      const unsigned int n_edges(edges.size());
      // Check that it is a new edge.
      ui_vector v_to_edge_pos(4,0);
      v_to_edge_pos[0] = vx_to_edge_pos.count(vertex_0[0]);
      v_to_edge_pos[1] = vy_to_edge_pos.count(vertex_0[1]);
      v_to_edge_pos[2] = vx_to_edge_pos.count(vertex_1[0]);
      v_to_edge_pos[3] = vy_to_edge_pos.count(vertex_1[1]);
      if (v_to_edge_pos[0]*v_to_edge_pos[1]*v_to_edge_pos[2]*v_to_edge_pos[3]>0)
      {
        pair<di_multimap::iterator,di_multimap::iterator> min_v;
        unsigned int pos(min_element(v_to_edge_pos.begin(),v_to_edge_pos.end())-
            v_to_edge_pos.begin());
        
        switch (pos)
        {
          case 0:
            {
              min_v = vx_to_edge_pos.equal_range(vertex_0[0]);
              break;
            }
          case 1:
            {
              min_v = vy_to_edge_pos.equal_range(vertex_0[1]);
              break;
            }
          case 2:
            {
              min_v = vx_to_edge_pos.equal_range(vertex_1[0]);
              break;
            }
          case 3:
            {
              min_v = vy_to_edge_pos.equal_range(vertex_1[1]);
              break;
            }
        }

        // Loop over the smallest set of edges
        di_multimap::iterator min_v_it(min_v.first);
        for (; min_v_it!=min_v.second; ++min_v_it)
          if (edges[min_v_it->second].Has_same_coord(vertex_0,vertex_1))
          {
            edges[min_v_it->second].Set_cell_index(1,i);
            cell_to_edge_gid[i].push_back(edges[min_v_it->second].Get_gid());
            if ((edges[min_v_it->second].Get_v0_x()==vertex_0[0]) &&
                (edges[min_v_it->second].Get_v0_y()==vertex_0[1]))
            {
              edges[min_v_it->second].Get_v0().push_back(vertex_0[2]);
              edges[min_v_it->second].Get_v1().push_back(vertex_1[2]);
            }
            else
            {
              edges[min_v_it->second].Get_v1().push_back(vertex_0[2]);
              edges[min_v_it->second].Get_v0().push_back(vertex_1[2]);
            }
            edge_exist = true;
            break;
          }
      }

      if (edge_exist==false)
      {
        EDGE_TYPE edge_type(Get_edge_type(vertex_0,vertex_1));
        edges.push_back(EDGE(edge_type,edge_gid,vertex_0,vertex_1));
        edges[n_edges].Set_cell_index(0,i);
        cell_to_edge_gid[i].push_back(edge_gid);
        // Store the coordinates of the vertices and the edge_gid associated
        // with the vertices
        vx_to_edge_pos.insert(di_pair(vertex_0[0],edge_gid));
        vy_to_edge_pos.insert(di_pair(vertex_0[1],edge_gid));
        vx_to_edge_pos.insert(di_pair(vertex_1[0],edge_gid));
        vy_to_edge_pos.insert(di_pair(vertex_1[1],edge_gid));
        ++edge_gid;
      }
    }
  }
}

unsigned int TRIANGULATION::Get_n_materials() const
{
  ui_set n_materials;
  const unsigned int mat_id_size(mat_id.size());
  for (unsigned int i=0; i<mat_id_size; ++i)
    n_materials.insert(mat_id[i]);

  return n_materials.size();
}

unsigned int TRIANGULATION::Get_n_sources() const
{
  ui_set n_sources;
  const unsigned int src_id_size(src_id.size());
  for (unsigned int i=0; i<src_id_size; ++i)
    n_sources.insert(src_id[i]);

  return n_sources.size();
}

EDGE_TYPE TRIANGULATION::Get_edge_type(d_vector const &vertex_0,
    d_vector const &vertex_1)
{
  EDGE_TYPE edge_type;
  if (vertex_0[1]==bottom_value && vertex_1[1]==bottom_value)
    edge_type = bottom_boundary;
  else
  {
    if (vertex_0[0]==right_value && vertex_1[0]==right_value)
      edge_type = right_boundary;
    else
    {
      if (vertex_0[1]==top_value && vertex_1[1]==top_value)
        edge_type = top_boundary;
      else
      {
        if (vertex_0[0]==left_value && vertex_1[0]==left_value)
          edge_type = left_boundary;
        else 
          edge_type = interior;
      }
    }
  }
  
  return edge_type;
}

void TRIANGULATION::Refine_mesh(ui_set const &cells_to_refine,
    ui_set const &adjacent_cells,
    map<unsigned int,vector<vector<d_vector> > > &edge_to_refine,
    vector<ui_vector> &projection)
{
  unsigned int new_n_cells(0);
  vector<vector<d_vector> > coarsest_vertices(n_cells);
  ui_vector new_mat_id,new_src_id,new_n_vertices,parents_map;
  ui_vector tmp_n_vertices(n_cells);
  vector<vector<d_vector> > new_grid;
  vector<vector<d_vector> > old_grid(grid);
  // Loop over the old cells
  for (unsigned int i=0; i<n_cells; ++i)
  {
    if (adjacent_cells.count(i)!=0)
    {
      // An adjacent cell is refined 
      unsigned int v(0);
      vector<d_vector> polygon;
      for (unsigned int j=0; j<n_vertices[i]; ++j)
      {
        // Add n existing vertex to the temporary cell
        d_vector vertex_0(grid[i][j]);
        d_vector vertex_1(grid[i][(j+1)%n_vertices[i]]);
        polygon.push_back(vertex_0);
        ++v;
        // If necessary create and add a new vertex to the temporary cell
        for (unsigned int k=0; k<edge_to_refine[i].size(); ++k)
        {
          if ((((vertex_0[0]==edge_to_refine[i][k][0][0]) && 
                (vertex_0[1]==edge_to_refine[i][k][0][1])) || 
               ((vertex_0[0]==edge_to_refine[i][k][1][0]) &&
                (vertex_0[1]==edge_to_refine[i][k][1][1]))) &&
              (((vertex_1[0]==edge_to_refine[i][k][0][0]) && 
                (vertex_1[1]==edge_to_refine[i][k][0][1])) || 
               ((vertex_1[0]==edge_to_refine[i][k][1][0]) &&
                (vertex_1[1]==edge_to_refine[i][k][1][1]))))
          {
            d_vector new_vertex(3,1.);
            if (j!=n_vertices[i]-1)
            {
              new_vertex[0] = 0.5*(grid[i][j][0]+grid[i][j+1][0]);
              new_vertex[1] = 0.5*(grid[i][j][1]+grid[i][j+1][1]);
            }
            else
            {
              new_vertex[0] = 0.5*(grid[i][j][0]+grid[i][0][0]);
              new_vertex[1] = 0.5*(grid[i][j][1]+grid[i][0][1]);
            }
            polygon.push_back(new_vertex);
            ++v;
          }
        }
      }
      grid[i] = polygon;
      tmp_n_vertices[i] = v;
    }
    else
      tmp_n_vertices[i] = n_vertices[i];
  }
        
  // Loop over the temporary cells
  for (unsigned int i=0; i<n_cells; ++i)
  {
    if (cells_to_refine.count(i)==0)
    {
      if (adjacent_cells.count(i)==0)
        new_n_vertices.push_back(n_vertices[i]);
      else
        new_n_vertices.push_back(tmp_n_vertices[i]);

      // If the cell is unchanged, we just add it to the new grid
      new_grid.push_back(grid[i]);
      new_mat_id.push_back(mat_id[i]);
      new_src_id.push_back(src_id[i]);
      parents_map.push_back(i);
      ++new_n_cells;
    }
    else
    {
      // The cell needs to be refined
      unsigned int n_coarser_vertices(0);
      double n_coarsest_vertices(0.);
      d_vector center(3,0.);
      vector<d_vector> coarse_vertices;
      // Extract the coarsest vertices
      for (unsigned int j=0; j<tmp_n_vertices[i]; ++j)
      {
        d_vector vertex(3,0.);
        vertex[0] = grid[i][j][0];
        vertex[1] = grid[i][j][1];
        coarsest_vertices[i].push_back(vertex);
        if (grid[i][j][2]==0.)
        {
          d_vector new_vertex(3,0.);
          coarse_vertices.push_back(new_vertex);
          coarse_vertices.push_back(vertex);
          center[0] += grid[i][j][0];
          center[1] += grid[i][j][1];
          n_coarsest_vertices += 1.;
        }
      }
      center[0] /= n_coarsest_vertices;
      center[1] /= n_coarsest_vertices;
      n_coarser_vertices = 2*static_cast<unsigned int>(n_coarsest_vertices);

      // Add the new coarser vertices
      for (unsigned int j=1; j<n_coarser_vertices; j+=2)
      {
        unsigned int jm2((j==1) ? n_coarser_vertices-1 : j-2);
        coarse_vertices[j-1][0] += 0.5*(coarse_vertices[j][0]+
            coarse_vertices[jm2][0]);
        coarse_vertices[j-1][1] += 0.5*(coarse_vertices[j][1]+
            coarse_vertices[jm2][1]);
        for (unsigned int k=0; k<grid[i].size(); ++k)
        {
          if ((grid[i][k][0]==coarse_vertices[j][0]) &&
              (grid[i][k][1]==coarse_vertices[j][1]))
          {
            unsigned int km1((k==0) ? (grid[i].size()-1) : (k-1));
            if ((grid[i][km1][0]==coarse_vertices[jm2][0]) &&
                (grid[i][km1][1]==coarse_vertices[jm2][1]))
              grid[i].insert(grid[i].begin()+k,coarse_vertices[j-1]);
          }
        }
      }

      unsigned int n_old_vertices(grid[i].size());
      for (unsigned int j=0; j<n_coarser_vertices; j+=2)
      {
        unsigned int pos(0);
        unsigned int jp2((j!=n_coarser_vertices-2) ? j+2 : 0);
        d_vector vertex(3,0.);
        vector<d_vector> polygon;
        // Add the first coarse vertex of the new cell
        vertex[0] = coarse_vertices[j][0];
        vertex[1] = coarse_vertices[j][1];
        polygon.push_back(vertex);
        // Find the previous coarsest vertex 
        while ((grid[i][pos][0]!=coarse_vertices[j][0]) || 
            (grid[i][pos][1]!=coarse_vertices[j][1]))
          pos = (pos+1)%n_old_vertices;
        // Move to the next vertex
        pos = (pos+1)%n_old_vertices;

        // Add the new vertices to the new cell
        while ((grid[i][pos][0]!=coarse_vertices[j+1][0]) ||
            (grid[i][pos][1]!=coarse_vertices[j+1][1]))
        {
          polygon.push_back(grid[i][pos]);
          pos = (pos+1)%n_old_vertices;
        }
                                            
        // Add a new coarse vertex
        vertex[0] = coarse_vertices[j+1][0];
        vertex[1] = coarse_vertices[j+1][1];
        polygon.push_back(vertex);
        pos = (pos+1)%n_old_vertices;

        // Add the new vertices to the new cell
        while ((grid[i][pos][0]!=coarse_vertices[jp2][0]) ||
            (grid[i][pos][1]!=coarse_vertices[jp2][1]))
        {
          polygon.push_back(grid[i][pos]);
          pos = (pos+1)%n_old_vertices;
        }

        // Add a new coarse vertex
        vertex[0] = coarse_vertices[jp2][0];
        vertex[1] = coarse_vertices[jp2][1];
        polygon.push_back(vertex);
        pos = (pos+1)%n_old_vertices;
        // Add the center of the cell
        polygon.push_back(center);
        new_grid.push_back(polygon);
        new_n_vertices.push_back(polygon.size());
        pos = (pos+1)%n_old_vertices;
        new_mat_id.push_back(mat_id[i]);
        new_src_id.push_back(src_id[i]);
        parents_map.push_back(i);
        ++new_n_cells;
      }
    }
  }
  grid = new_grid;
  projection = Project(old_grid,grid,parents_map,new_n_cells);
  n_vertices = new_n_vertices;
  mat_id = new_mat_id;
  src_id = new_src_id;
  n_cells = new_n_cells;
  edges.clear();
  cell_to_edge_gid.clear();
}

vector<ui_vector> TRIANGULATION::Project(vector<vector<d_vector> > const &old_grid,
    vector<vector<d_vector> > const &new_grid,ui_vector const &parents_map,
    const unsigned int n_cells)
{
  unsigned int dof(0);
  vector<ui_vector> projection;
  // Loop over the cells of the new grids
  for (unsigned int i=0; i<n_cells; ++i)
  {
    const unsigned int n_vertices(new_grid[i].size());
    const unsigned int parents_i(parents_map[i]);
    const unsigned int n_old_vertices(old_grid[parents_i].size());
    // Loop over the vertices of the cell
    for (unsigned int j=0; j<n_vertices; ++j)
    {
      bool existing_vertex(false);
      ui_vector vertex_proj;

      // Check if the current vertex existed in the old grid
      for (unsigned int k=0; k<n_old_vertices; ++k)
      {
        if ((grid[i][j][0]==old_grid[parents_i][k][0]) &&
            (grid[i][j][1]==old_grid[parents_i][k][1]))
        {
          vertex_proj.push_back(dof+k);
          existing_vertex = true;
          break;
        }
      }

      // If the vertex did not exist, the flux will have to be interpolated
      if (existing_vertex==false)
      {
        bool center_vertex(true);
        for (unsigned int k=0; k<n_old_vertices; ++k)
        {
          const unsigned int kp1((k+1)%n_old_vertices);
          double new_x = 0.5*(old_grid[parents_i][k][0]+old_grid[parents_i][kp1][0]);
          double new_y = 0.5*(old_grid[parents_i][k][1]+old_grid[parents_i][kp1][1]);
          if ((grid[i][j][0]==new_x) && (grid[i][j][1]==new_y))
          {
            vertex_proj.push_back(dof+k);
            vertex_proj.push_back(dof+kp1);
            center_vertex = false;
            break;
          }
        }
        if (center_vertex==true)
          for (unsigned int k=0; k<n_old_vertices; ++k)
            vertex_proj.push_back(dof+k);
      }
      projection.push_back(vertex_proj);
    }
    if (i!=n_cells-1)
      if (parents_i!=parents_map[i+1])
        dof += n_old_vertices;
  }

  return projection;
}
