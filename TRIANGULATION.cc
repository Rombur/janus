#include "TRIANGULATION.hh"

void TRIANGULATION::Read_geometry()
{
  // Open the file to read it.
  ifstream geometry_file(geometry_inputfile,ios::in);

  // Read the cell type: quadrilateral or polygon
  geometry_file>>cell_type;

  const unsigned int dim(2);
  if (cell_type==rectangle)
  {
    // Read the number of cells in x and y
    geometry_file>>n_x>>n_y;
    n_cells = n_x*n_y;
    n_vertices.resize(n_cells,4);
    // Read the abscissae
    const unsigned int n_x_vertices(n_x+1);
    for (unsigned int i=0; i<n_x_vertices; ++i)
      geometry_file>>x[i];
    // Read the ordinates
    const unsigned int n_y_vertices(n_y+1);
    for (unsigned int i=0; i<n_y_vertices; ++i)
      geometry_file>>y[i];
    // Read the material ids
    for (unsigned int i=0; i<n_cells; ++i)
      geometry_file>>mat_id[i];
    // Read the source ids
    for (unsigned int i=0; i<n_cells; ++i)
      geometry_file>>src_id[i];
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
        polygon[vert] = coord;
      }
      grid[cell] = polygon;
      // Read the material id
      geometry_file>>mat_id[cell];
      // Read the source id
      geometry_file>>src_id[cell];
    }
  }
}

