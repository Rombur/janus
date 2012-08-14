#include "TRIANGULATION.hh"

TRIANGULATION::TRIANGULATION(string* geometry_inputfile) : 
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
  assert(geometry_file);

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
    unsigned int n_x(0),n_y(0);
    // Read the number of cells in x and y
    geometry_file>>n_x>>n_y;
    n_cells = n_x*n_y;
    mat_id.resize(n_cells);
    src_id.resize(n_cells);
    // Fill n_vertices with 4 since all the cells are rectangular
    n_vertices.resize(n_cells,4);
    // Read the abscissae
    const unsigned int n_x_vertices(n_x+1);
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
    const unsigned int n_y_vertices(n_y+1);
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
      vector<d_vector> rectangle(4,d_vector(dim,0.));
      if (i%n_x==0 && i!=0)
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
        d_vector coord(dim,0.);
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
      d_vector vertex_0(dim,0.);
      d_vector vertex_1(dim,0.);
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
