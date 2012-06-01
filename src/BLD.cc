#include "BLD.hh"

BLD::BLD(d_vector const &cell_x,d_vector const &cell_y) :
  FINITE_ELEMENT(cell_x,cell_y)
{
  dof_per_cell = 4;
  delta_x = x[1]-x[0];
  delta_y = y[3]-y[0];
  Teuchos::SerialDenseMatrix<int,double> zeros_4(4,4);
  downwind.resize(4,zeros_4);
  upwind.resize(4,zeros_4);
  edge_deln_matrix.resize(4,vector<Teuchos::SerialDenseMatrix<int,double> > (2,
        zeros_4));
  coupling_edge_deln_matrix.resize(4,
      vector<Teuchos::SerialDenseMatrix<int,double> > (2,zeros_4));
  mass_matrix = zeros_4;
  x_grad_matrix = zeros_4;
  y_grad_matrix = zeros_4;
  stiffness_matrix = zeros_4;
  x_grad_edge_matrix.resize(4,zeros_4);
  y_grad_edge_matrix.resize(4,zeros_4);
}

void BLD::Build_fe_1d()
{
  double h_ratio(delta_x/(6.*delta_y));
  double v_ratio(delta_y/(6.*delta_x));

  vector<d_vector> horizontal_edge_mass_matrix(2,d_vector(2,0.));
  vector<d_vector> vertical_edge_mass_matrix(2,d_vector(2,0.));

  horizontal_edge_mass_matrix[0][0] = delta_x/3.;
  horizontal_edge_mass_matrix[0][1] = delta_x/6.;
  horizontal_edge_mass_matrix[1][0] = delta_x/6.;
  horizontal_edge_mass_matrix[1][1] = delta_x/3.;

  vertical_edge_mass_matrix[0][0] = delta_y/3.;
  vertical_edge_mass_matrix[0][1] = delta_y/6.;
  vertical_edge_mass_matrix[1][0] = delta_y/6.;
  vertical_edge_mass_matrix[1][1] = delta_y/3.;

  // Build the downwind matrices
  // Matrix associated to the bottom edge
  downwind[0](0,0) = horizontal_edge_mass_matrix[0][0];
  downwind[0](0,1) = horizontal_edge_mass_matrix[0][1];
  downwind[0](1,0) = horizontal_edge_mass_matrix[1][0];
  downwind[0](1,1) = horizontal_edge_mass_matrix[1][1];
  // Matrix associated to the right edge
  downwind[1](1,1) = vertical_edge_mass_matrix[0][0];
  downwind[1](1,2) = vertical_edge_mass_matrix[0][1];
  downwind[1](2,1) = vertical_edge_mass_matrix[1][0];
  downwind[1](2,2) = vertical_edge_mass_matrix[1][1];
  // Matrix associated to the top edge
  downwind[2](2,2) = horizontal_edge_mass_matrix[0][0];
  downwind[2](2,3) = horizontal_edge_mass_matrix[0][1];
  downwind[2](3,2) = horizontal_edge_mass_matrix[1][0];
  downwind[2](3,3) = horizontal_edge_mass_matrix[1][1];
  // Matrix associated to the left edge
  downwind[3](0,0) = vertical_edge_mass_matrix[0][0];
  downwind[3](0,3) = vertical_edge_mass_matrix[0][1];
  downwind[3](3,0) = vertical_edge_mass_matrix[1][0];
  downwind[3](3,3) = vertical_edge_mass_matrix[1][1];

  // Build the edge_deln_matrices
  // Matrix associated to the bottom edge
  edge_deln_matrix[0][1](0,0) = -2.*h_ratio;
  edge_deln_matrix[0][1](0,1) = -h_ratio;
  edge_deln_matrix[0][1](1,0) = -h_ratio;
  edge_deln_matrix[0][1](1,1) = -2.*h_ratio;
  edge_deln_matrix[0][1](2,0) = h_ratio;
  edge_deln_matrix[0][1](2,1) = 2.*h_ratio;
  edge_deln_matrix[0][1](3,0) = 2.*h_ratio;
  edge_deln_matrix[0][1](3,1) = h_ratio;
  // Matrix associated to the right edge
  edge_deln_matrix[1][0](0,1) = -2.*v_ratio;
  edge_deln_matrix[1][0](0,2) = -v_ratio;
  edge_deln_matrix[1][0](1,1) = 2.*v_ratio;
  edge_deln_matrix[1][0](1,2) = v_ratio;
  edge_deln_matrix[1][0](2,1) = v_ratio;
  edge_deln_matrix[1][0](2,2) = 2.*v_ratio;
  edge_deln_matrix[1][0](3,1) = -v_ratio;
  edge_deln_matrix[1][0](3,2) = -2.*v_ratio;
  // Matrix associated to the top edge
  edge_deln_matrix[2][1](0,2) = -h_ratio;
  edge_deln_matrix[2][1](0,3) = -2.*h_ratio;
  edge_deln_matrix[2][1](1,2) = -2.*h_ratio;
  edge_deln_matrix[2][1](1,3) = -h_ratio;
  edge_deln_matrix[2][1](2,2) = 2.*h_ratio;
  edge_deln_matrix[2][1](2,3) = h_ratio;
  edge_deln_matrix[2][1](3,2) = h_ratio;
  edge_deln_matrix[2][1](3,3) = 2.*h_ratio;
  // Matrix associated to the left edge
  edge_deln_matrix[3][0](0,0) = -2.*v_ratio;
  edge_deln_matrix[3][0](0,3) = -v_ratio;
  edge_deln_matrix[3][0](1,0) = 2.*v_ratio;
  edge_deln_matrix[3][0](1,3) = v_ratio;
  edge_deln_matrix[3][0](2,0) = v_ratio;
  edge_deln_matrix[3][0](2,3) = 2.*v_ratio;
  edge_deln_matrix[3][0](3,0) = -v_ratio;
  edge_deln_matrix[3][0](3,3) = -2.*v_ratio;

  // Build the grad_edge_matrices
  // Matrix associated to the bottom edge
  x_grad_edge_matrix[0](0,0) = -1./delta_x;
  x_grad_edge_matrix[0](0,1) = 1./delta_x;
  x_grad_edge_matrix[0](1,0) = -1./delta_x;
  x_grad_edge_matrix[0](1,1) = 1./delta_x;
  y_grad_edge_matrix[0](0,0) = -1./delta_y;
  y_grad_edge_matrix[0](0,3) = 1./delta_y;
  y_grad_edge_matrix[0](1,1) = -1./delta_y;
  y_grad_edge_matrix[0](1,2) = -1./delta_y;
  // Matrix associated to the right edge
  x_grad_edge_matrix[1](1,0) = -1./delta_x;
  x_grad_edge_matrix[1](1,1) = 1./delta_x;
  x_grad_edge_matrix[1](2,2) = 1./delta_x;
  x_grad_edge_matrix[1](2,3) = -1./delta_x;
  y_grad_edge_matrix[1](1,1) = -1./delta_y;
  y_grad_edge_matrix[1](1,2) = 1./delta_y;
  y_grad_edge_matrix[1](2,1) = -1./delta_y;
  y_grad_edge_matrix[1](2,2) = 1./delta_y;
  // Matrix associated to the top edge
  x_grad_edge_matrix[2](2,2) = 1./delta_x;
  x_grad_edge_matrix[2](2,3) = -1./delta_x;
  x_grad_edge_matrix[2](3,2) = 1./delta_x;
  x_grad_edge_matrix[2](3,3) = -1./delta_x;
  y_grad_edge_matrix[2](2,1) = -1./delta_y;
  y_grad_edge_matrix[2](2,2) = 1./delta_y;
  y_grad_edge_matrix[2](3,0) = -1./delta_y;
  y_grad_edge_matrix[2](3,3) = 1./delta_y;
  // Matrix associated to the right edge
  x_grad_edge_matrix[3](0,0) = -1./delta_x;
  x_grad_edge_matrix[3](0,1) = 1./delta_x;
  x_grad_edge_matrix[3](3,2) = 1./delta_x;
  x_grad_edge_matrix[3](3,3) = -1./delta_x;
  y_grad_edge_matrix[3](0,0) = -1./delta_y;
  y_grad_edge_matrix[3](0,3) = 1./delta_y;
  y_grad_edge_matrix[3](3,0) = -1./delta_y;
  y_grad_edge_matrix[3](3,3) = 1./delta_y;
}

void BLD::Build_upwind_matrices(CELL* cell,vector<CELL*> const &mesh)
{
  double h_ratio(delta_x/(6.*delta_y));
  double v_ratio(delta_y/(6.*delta_x));

  vector<d_vector> coupling_horizontal_edge_mass_matrix(2,d_vector(2,0.));
  vector<d_vector> coupling_vertical_edge_mass_matrix(2,d_vector(2,0.));

  coupling_horizontal_edge_mass_matrix[0][0] = delta_x/6.;
  coupling_horizontal_edge_mass_matrix[0][1] = delta_x/3.;
  coupling_horizontal_edge_mass_matrix[1][0] = delta_x/3.;
  coupling_horizontal_edge_mass_matrix[1][1] = delta_x/6.;

  coupling_vertical_edge_mass_matrix[0][0] = delta_y/3.;
  coupling_vertical_edge_mass_matrix[0][1] = delta_y/6.;
  coupling_vertical_edge_mass_matrix[1][0] = delta_y/6.;
  coupling_vertical_edge_mass_matrix[1][1] = delta_y/3.;

  // Build the upwind matrices
  // Matrix associated to the bottom edge
  upwind[0](0,2) = coupling_horizontal_edge_mass_matrix[0][0];
  upwind[0](0,3) = coupling_horizontal_edge_mass_matrix[0][1];
  upwind[0](1,2) = coupling_horizontal_edge_mass_matrix[1][0];
  upwind[0](1,3) = coupling_horizontal_edge_mass_matrix[1][1];
  // Matrix associated to the right edge
  upwind[1](1,0) = coupling_vertical_edge_mass_matrix[0][0];
  upwind[1](1,3) = coupling_vertical_edge_mass_matrix[0][1];
  upwind[1](2,0) = coupling_vertical_edge_mass_matrix[1][0];
  upwind[1](2,3) = coupling_vertical_edge_mass_matrix[1][1];
  // Matrix associated to the top edge
  upwind[2](2,0) = coupling_horizontal_edge_mass_matrix[0][0];
  upwind[2](2,1) = coupling_horizontal_edge_mass_matrix[0][1];
  upwind[2](3,0) = coupling_horizontal_edge_mass_matrix[1][0];
  upwind[2](3,1) = coupling_horizontal_edge_mass_matrix[1][1];
  // Matrix associated to the left edge
  upwind[3](0,1) = coupling_vertical_edge_mass_matrix[0][0];
  upwind[3](0,2) = coupling_vertical_edge_mass_matrix[0][1];
  upwind[3](3,1) = coupling_vertical_edge_mass_matrix[1][0];
  upwind[3](3,2) = coupling_vertical_edge_mass_matrix[1][1];

  // Build the edge_deln_matrices
  // Matrix associated to the bottom edge
  coupling_edge_deln_matrix[0][1](0,2) = -h_ratio;
  coupling_edge_deln_matrix[0][1](0,3) = -2.*h_ratio;
  coupling_edge_deln_matrix[0][1](1,2) = -2.*h_ratio;
  coupling_edge_deln_matrix[0][1](1,3) = -h_ratio;
  coupling_edge_deln_matrix[0][1](2,2) = 2.*h_ratio;
  coupling_edge_deln_matrix[0][1](2,3) = h_ratio;
  coupling_edge_deln_matrix[0][1](3,2) = h_ratio;
  coupling_edge_deln_matrix[0][1](3,3) = 2.*h_ratio;
  // Matrix associated to the right edge
  coupling_edge_deln_matrix[1][0](0,0) = -2.*v_ratio;
  coupling_edge_deln_matrix[1][0](0,3) = -v_ratio;
  coupling_edge_deln_matrix[1][0](1,0) = 2.*v_ratio;
  coupling_edge_deln_matrix[1][0](1,3) = v_ratio;
  coupling_edge_deln_matrix[1][0](2,0) = v_ratio;
  coupling_edge_deln_matrix[1][0](2,3) = 2.*v_ratio;
  coupling_edge_deln_matrix[1][0](3,0) = -v_ratio;
  coupling_edge_deln_matrix[1][0](3,3) = -2.*v_ratio;
  // Matrix associated to the top edge
  coupling_edge_deln_matrix[2][1](0,0) = -2.*h_ratio;
  coupling_edge_deln_matrix[2][1](0,1) = -h_ratio;
  coupling_edge_deln_matrix[2][1](1,0) = -h_ratio;
  coupling_edge_deln_matrix[2][1](1,1) = -2.*h_ratio;
  coupling_edge_deln_matrix[2][1](2,0) = h_ratio;
  coupling_edge_deln_matrix[2][1](2,1) = 2.*h_ratio;
  coupling_edge_deln_matrix[2][1](3,0) = 2.*h_ratio;
  coupling_edge_deln_matrix[2][1](3,1) = h_ratio;  
  // Matrix associated to the left edge
  coupling_edge_deln_matrix[3][0](0,1) = -2.*v_ratio;
  coupling_edge_deln_matrix[3][0](0,2) = -v_ratio;
  coupling_edge_deln_matrix[3][0](1,1) = 2.*v_ratio;
  coupling_edge_deln_matrix[3][0](1,2) = v_ratio;
  coupling_edge_deln_matrix[3][0](2,1) = v_ratio;
  coupling_edge_deln_matrix[3][0](2,2) = 2.*v_ratio;
  coupling_edge_deln_matrix[3][0](3,1) = -v_ratio;
  coupling_edge_deln_matrix[3][0](3,2) = -2.*v_ratio;
}

void BLD::Build_fe_2d()
{
  const double mass_ratio(delta_x*delta_y/36.);
  const double x_grad_ratio(delta_y/12.);
  const double y_grad_ratio(delta_x/12.);
  const double x_stiffness_ratio(delta_y/(6.*delta_x));
  const double y_stiffness_ratio(delta_x/(6.*delta_y));

  // Build the mass matrix
  mass_matrix(0,0) = 4.*mass_ratio;
  mass_matrix(0,1) = 2.*mass_ratio;
  mass_matrix(0,2) = mass_ratio;
  mass_matrix(0,3) = 2.*mass_ratio;
  mass_matrix(1,0) = 2.*mass_ratio;
  mass_matrix(1,1) = 4.*mass_ratio;
  mass_matrix(1,2) = 2.*mass_ratio;
  mass_matrix(1,3) = mass_ratio;
  mass_matrix(2,0) = mass_ratio;
  mass_matrix(2,1) = 2.*mass_ratio;
  mass_matrix(2,2) = 4.*mass_ratio;
  mass_matrix(2,3) = 2.*mass_ratio;
  mass_matrix(3,0) = 2.*mass_ratio;
  mass_matrix(3,1) = mass_ratio;
  mass_matrix(3,2) = 2.*mass_ratio;
  mass_matrix(3,3) = 4.*mass_ratio;

  // Build the x_grad matrix
  x_grad_matrix(0,0) = -2.*x_grad_ratio;
  x_grad_matrix(0,1) = -2.*x_grad_ratio;
  x_grad_matrix(0,2) = -x_grad_ratio;
  x_grad_matrix(0,3) = -x_grad_ratio;
  x_grad_matrix(1,0) = 2.*x_grad_ratio;
  x_grad_matrix(1,1) = 2.*x_grad_ratio;
  x_grad_matrix(1,2) = x_grad_ratio;
  x_grad_matrix(1,3) = x_grad_ratio;
  x_grad_matrix(2,0) = x_grad_ratio;
  x_grad_matrix(2,1) = x_grad_ratio;
  x_grad_matrix(2,2) = 2.*x_grad_ratio;
  x_grad_matrix(2,3) = 2.*x_grad_ratio;
  x_grad_matrix(3,0) = -x_grad_ratio;
  x_grad_matrix(3,1) = -x_grad_ratio;
  x_grad_matrix(3,2) = -2.*x_grad_ratio;
  x_grad_matrix(3,3) = -2.*x_grad_ratio;

  // Build the y_grad matrix
  y_grad_matrix(0,0) = -2.*y_grad_ratio;
  y_grad_matrix(0,1) = -y_grad_ratio;
  y_grad_matrix(0,2) = -y_grad_ratio;
  y_grad_matrix(0,3) = -2.*y_grad_ratio;
  y_grad_matrix(1,0) = -y_grad_ratio;
  y_grad_matrix(1,1) = -2.*y_grad_ratio;
  y_grad_matrix(1,2) = -2.*y_grad_ratio;
  y_grad_matrix(1,3) = -y_grad_ratio;
  y_grad_matrix(2,0) = y_grad_ratio;
  y_grad_matrix(2,1) = 2.*y_grad_ratio;
  y_grad_matrix(2,2) = 2.*y_grad_ratio;
  y_grad_matrix(2,3) = y_grad_ratio;
  y_grad_matrix(3,0) = 2.*y_grad_ratio;
  y_grad_matrix(3,1) = y_grad_ratio;
  y_grad_matrix(3,2) = y_grad_ratio;
  y_grad_matrix(3,3) = 2.*y_grad_ratio;

  // Build the stiffness matrix
  stiffness_matrix(0,0) =  2.*x_stiffness_ratio+2.*y_stiffness_ratio;
  stiffness_matrix(0,1) =  -2.*x_stiffness_ratio+y_stiffness_ratio;
  stiffness_matrix(0,2) =  -x_stiffness_ratio-y_stiffness_ratio;
  stiffness_matrix(0,3) =  x_stiffness_ratio-2.*y_stiffness_ratio;
  stiffness_matrix(1,0) =  -2.*x_stiffness_ratio+y_stiffness_ratio;
  stiffness_matrix(1,1) =  2.*x_stiffness_ratio+2.*y_stiffness_ratio;
  stiffness_matrix(1,2) =  x_stiffness_ratio-2.*y_stiffness_ratio;
  stiffness_matrix(1,3) =  -x_stiffness_ratio-y_stiffness_ratio;
  stiffness_matrix(2,0) =  -x_stiffness_ratio-y_stiffness_ratio;
  stiffness_matrix(2,1) =  x_stiffness_ratio-2.*y_stiffness_ratio;
  stiffness_matrix(2,2) =  2.*x_stiffness_ratio+2.*y_stiffness_ratio;
  stiffness_matrix(2,3) =  -2.*x_stiffness_ratio+y_stiffness_ratio;
  stiffness_matrix(3,0) =  x_stiffness_ratio-2.*y_stiffness_ratio;
  stiffness_matrix(3,1) =  -x_stiffness_ratio-y_stiffness_ratio;
  stiffness_matrix(3,2) =  -2.*x_stiffness_ratio+y_stiffness_ratio;
  stiffness_matrix(3,3) =  2.*x_stiffness_ratio+2.*y_stiffness_ratio;
}
