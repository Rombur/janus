#include "PWLD.hh"

PWLD::PWLD(d_vector const &cell_x,d_vector const &cell_y) :
  FINITE_ELEMENT(cell_x,cell_y)
{
  const unsigned int dim(2);
  dof_per_cell = x.size();
  x_c = 0.;
  y_c = 0.;
  for (unsigned int i=0; i<dof_per_cell; i++)
  {
    x_c += x[i];
    y_c += y[i];
  }
  x_c /= dof_per_cell;
  y_c /= dof_per_cell;
  Teuchos::SerialDenseMatrix<int,double> zeros(dof_per_cell,dof_per_cell);
  downwind.resize(dof_per_cell,zeros);
  upwind.resize(dof_per_cell);
  edge_deln_matrix.resize(dof_per_cell,
      vector<Teuchos::SerialDenseMatrix<int,double> > (dim,zeros));
  coupling_edge_deln_matrix.resize(dof_per_cell,
      vector<Teuchos::SerialDenseMatrix<int,double> > (dim));
  mass_matrix = zeros;
  x_grad_matrix = zeros;
  y_grad_matrix = zeros;
  stiffness_matrix = zeros;
}

void PWLD::Build_fe_1d()
{    
  const unsigned int dim(2);
  vector<d_vector> array_0(2,d_vector(2,0.));
  vector<d_vector> array_1(3,d_vector(3,0.));
  vector<d_vector> array_2(3,d_vector(3,0.));
  array_0[0][0] = 2.;
  array_0[0][1] = 1.;
  array_0[1][0] = 1.;
  array_0[1][1] = 2.;
  array_1[0][0] = -1.;
  array_1[0][1] = -1.;
  array_1[1][0] = 1.;
  array_1[1][1] = 1.;
  array_2[0][0] = -1.;
  array_2[0][1] = -1.;
  array_2[2][0] = 1.;
  array_2[2][1] = 1.;
  deln_matrix.resize(dof_per_cell,vector<vector<d_vector> > (2));
  for (unsigned int i=0; i<dof_per_cell; ++i)
  {
    double x0(x[i]);
    double y0(y[i]);
    double x1,y1,x1_x0,x2_x0,y1_y0,y2_y0,length,jacobian;
    vector<d_vector> tmp_array_0(2,d_vector(2,0.));
    vector<d_vector> tmp_array_1(3,d_vector(3,0.));
    vector<d_vector> tmp_array_2(3,d_vector(3,0.));
    if (i==dof_per_cell-1)
    {
      x1 = x[0];
      y1 = y[0];
    }
    else
    {
      x1 = x[i+1];
      y1 = y[i+1];
    }
    x1_x0 = x1-x0;
    x2_x0 = x_c-x0;
    y1_y0 = y1-y0;
    y2_y0 = y_c-y0;
    length = sqrt(pow(x1_x0,2)+pow(y1_y0,2));
    jacobian = fabs(x1_x0*y2_y0-y1_y0*x2_x0);
    for (unsigned int j=0; j<2; ++j)
      for (unsigned int k=0; k<2; ++k)
        tmp_array_0[j][k] = length/6.*array_0[j][k];
    for (unsigned int j=0; j<3; ++j)
      for (unsigned int k=0; k<3; ++k)
      {
        tmp_array_1[j][k] = length/(2.*jacobian)*(y2_y0*array_1[j][k]-
            y1_y0*array_2[j][k]);
        tmp_array_2[j][k] = length/(2.*jacobian)*(-x2_x0*array_1[j][k]+
            x1_x0*array_2[j][k]);
      }
    edge_mass_matrix.push_back(tmp_array_0);
    deln_matrix[i][0] = tmp_array_1;
    deln_matrix[i][1] = tmp_array_2;
  }

  for (unsigned int index_0=0; index_0<dof_per_cell; ++index_0)
  {
    unsigned int index_1((index_0+1)%dof_per_cell);
    // Build the downwind matrices
    downwind[index_0](index_0,index_0) = edge_mass_matrix[index_0][0][0];
    downwind[index_0](index_0,index_1) = edge_mass_matrix[index_0][0][1];
    downwind[index_0](index_1,index_0) = edge_mass_matrix[index_0][1][0];
    downwind[index_0](index_1,index_1) = edge_mass_matrix[index_0][1][1];
    // Build the edge deln matrices
    for (unsigned int i=0; i<dim; ++i)
    {
      edge_deln_matrix[index_0][i](index_0,index_0) = 
        deln_matrix[index_0][i][0][0];
      edge_deln_matrix[index_0][i](index_0,index_1) = 
        deln_matrix[index_0][i][0][1];
      edge_deln_matrix[index_0][i](index_1,index_0) = 
        deln_matrix[index_0][i][1][0];
      edge_deln_matrix[index_0][i](index_1,index_1) = 
        deln_matrix[index_0][i][1][1];
      for (unsigned int j=0; j<dof_per_cell; ++j)
      {
        edge_deln_matrix[index_0][i](j,index_0) += 1./dof_per_cell*
          deln_matrix[index_0][i][2][0];
        edge_deln_matrix[index_0][i](j,index_1) += 1./dof_per_cell*
          deln_matrix[index_0][i][2][1];
      }
    }
  }
}

void PWLD::Build_upwind_matrices(CELL* cell,vector<CELL*> const &mesh)
{
  const unsigned int dim(2);
  vector<EDGE*>::iterator cell_edge(cell->Get_cell_edges_begin());
  vector<EDGE*>::const_iterator cell_edges_end(cell->Get_cell_edges_end());

  for (unsigned int i=0; cell_edge<cell_edges_end; ++cell_edge,++i)
  {
    // The edge is on the boundary, there is no upwinding
    if ((*cell_edge)->Is_interior()==true)
    {
      unsigned int edge_lid(0);
      unsigned int upwind_edge_lid(0);
      unsigned int upwind_cell_id(0);
      unsigned int upwind_dof_per_cell(0);
      if (cell->Get_id()==(*cell_edge)->Get_cell_index(0))
      {
        edge_lid = (*cell_edge)->Get_lid(0);
        upwind_edge_lid = (*cell_edge)->Get_lid(1);
        upwind_cell_id = (*cell_edge)->Get_cell_index(1);
      }
      else
      {
        edge_lid = (*cell_edge)->Get_lid(1);
        upwind_edge_lid = (*cell_edge)->Get_lid(0);
        upwind_cell_id = (*cell_edge)->Get_cell_index(0);
      }
      // Read right to left: const pointer to a const FINITE_ELEMENT (pointer
      // and the value cannot be changed)
      FINITE_ELEMENT const* const upwind_fe(mesh[upwind_cell_id]->Get_fe());
      upwind_dof_per_cell = upwind_fe->Get_dof_per_cell();
      // Build the upwind matrix associated to the edge
      upwind[i].shape(dof_per_cell,upwind_dof_per_cell);
      upwind[i](edge_lid,(upwind_edge_lid+1)%upwind_dof_per_cell) +=
        edge_mass_matrix[i][0][0];
      upwind[i](edge_lid,upwind_edge_lid) += edge_mass_matrix[i][0][1];
      upwind[i]((edge_lid+1)%dof_per_cell,upwind_edge_lid) +=
        edge_mass_matrix[i][1][1];
      upwind[i]((edge_lid+1)%dof_per_cell,
        (upwind_edge_lid+1)%upwind_dof_per_cell) += 
          edge_mass_matrix[i][1][0];
      // Build the coupling_edge_deln_matrices associated to the edge
      for (unsigned int j=0; j<dim; ++j)
      {
        coupling_edge_deln_matrix[i][j].shape(dof_per_cell,upwind_dof_per_cell);
        coupling_edge_deln_matrix[i][j](edge_lid,upwind_edge_lid) +=
          deln_matrix[i][j][0][0];
        coupling_edge_deln_matrix[i][j](edge_lid,(upwind_edge_lid+1)%
          upwind_dof_per_cell) += deln_matrix[i][j][0][1];
        coupling_edge_deln_matrix[i][j]((edge_lid+1)%dof_per_cell,upwind_edge_lid) +=
          deln_matrix[i][j][1][0];
        coupling_edge_deln_matrix[i][j]((edge_lid+1)%dof_per_cell,
          (upwind_edge_lid+1)%upwind_dof_per_cell) += deln_matrix[i][j][1][1];
      }
    }
  }

  // Delete the temporary vector.
  edge_mass_matrix.clear();
  deln_matrix.clear();
}

void PWLD::Build_fe_2d()
{
  const double ratio(1./dof_per_cell);
  const double square_ratio(1./(dof_per_cell*dof_per_cell));
  // Build the matrices by looping over the sub-triangle
  for (unsigned int i=0; i<dof_per_cell; ++i)
  {
    unsigned int j((i+1)%dof_per_cell);
    double x0(x[i]);
    double x1(x[j]);
    double y0(y[i]);
    double y1(y[j]);
    double x1_x0(x1-x0);
    double x2_x1(x_c-x1);
    double x2_x0(x_c-x0);
    double y1_y0(y1-y0);
    double y2_y1(y_c-y1);
    double y2_y0(y_c-y0);
    double a_00(0.5*(pow(y2_y1,2)+pow(x2_x1,2)));
    double a_01(-0.5*(y2_y0*y2_y1+x2_x0*x2_x1));
    double a_02(0.5*(y1_y0*y2_y1+x1_x0*x2_x1));
    double a_11(0.5*(pow(y2_y0,2)+pow(x2_x0,2)));
    double a_12(-0.5*(y2_y0*y1_y0+x2_x0*x1_x0));
    double a_22(0.5*(pow(y1_y0,2)+pow(x1_x0,2)));
    double jacobian(fabs(x1_x0*y2_y0-y1_y0*x2_x0));
    vector<d_vector> mass(3,d_vector(3,0.));
    vector<d_vector> x_grad(3,d_vector(3,0.));
    vector<d_vector> y_grad(3,d_vector(3,0.));
    vector<d_vector> stiffness(3,d_vector(3,0.));

    mass[0][0] = jacobian/12.;
    mass[0][1] = jacobian/24.;
    mass[0][2] = jacobian/24.;
    mass[1][0] = jacobian/24.;
    mass[1][1] = jacobian/12.;
    mass[1][2] = jacobian/24.;
    mass[2][0] = jacobian/24.;
    mass[2][1] = jacobian/24.;
    mass[2][2] = jacobian/12.;

    x_grad[0][0] = -y2_y1/6.;
    x_grad[0][1] = -y2_y1/6.;
    x_grad[0][2] = -y2_y1/6.;
    x_grad[1][0] = y2_y0/6.;
    x_grad[1][1] = y2_y0/6.;
    x_grad[1][2] = y2_y0/6.;
    x_grad[2][0] = -y1_y0/6.;
    x_grad[2][1] = -y1_y0/6.;
    x_grad[2][2] = -y1_y0/6.;

    y_grad[0][0] = x2_x1/6.; 
    y_grad[0][1] = x2_x1/6.; 
    y_grad[0][2] = x2_x1/6.; 
    y_grad[1][0] = -x2_x0/6.;
    y_grad[1][1] = -x2_x0/6.;
    y_grad[1][2] = -x2_x0/6.;
    y_grad[2][0] = x1_x0/6.;
    y_grad[2][1] = x1_x0/6.;
    y_grad[2][2] = x1_x0/6.;

    stiffness[0][0] = a_00/jacobian;
    stiffness[0][1] = a_01/jacobian;
    stiffness[0][2] = a_02/jacobian;
    stiffness[1][0] = a_01/jacobian;
    stiffness[1][1] = a_11/jacobian;
    stiffness[1][2] = a_12/jacobian;
    stiffness[2][0] = a_02/jacobian;
    stiffness[2][1] = a_12/jacobian;
    stiffness[2][2] = a_22/jacobian;

    mass_matrix(i,i) += mass[0][0];
    mass_matrix(i,j) += mass[0][1];
    mass_matrix(j,i) += mass[1][0];
    mass_matrix(j,j) += mass[1][1];

    x_grad_matrix(i,i) += x_grad[0][0];
    x_grad_matrix(i,j) += x_grad[0][1];
    x_grad_matrix(j,i) += x_grad[1][0];
    x_grad_matrix(j,j) += x_grad[1][1];

    y_grad_matrix(i,i) += y_grad[0][0];
    y_grad_matrix(i,j) += y_grad[0][1];
    y_grad_matrix(j,i) += y_grad[1][0];
    y_grad_matrix(j,j) += y_grad[1][1];

    stiffness_matrix(i,i) += stiffness[0][0];
    stiffness_matrix(i,j) += stiffness[0][1];
    stiffness_matrix(j,i) += stiffness[1][0];
    stiffness_matrix(j,j) += stiffness[1][1];

    for (unsigned int k=0; k<dof_per_cell; ++k)
    {
      mass_matrix(i,k) += ratio*mass[0][2];
      mass_matrix(j,k) += ratio*mass[1][2];
      mass_matrix(k,i) += ratio*mass[2][0];
      mass_matrix(k,j) += ratio*mass[2][1];

      x_grad_matrix(i,k) += ratio*x_grad[0][2];
      x_grad_matrix(j,k) += ratio*x_grad[1][2];
      x_grad_matrix(k,i) += ratio*x_grad[2][0];
      x_grad_matrix(k,j) += ratio*x_grad[2][1];

      y_grad_matrix(i,k) += ratio*y_grad[0][2];
      y_grad_matrix(j,k) += ratio*y_grad[1][2];
      y_grad_matrix(k,i) += ratio*y_grad[2][0];
      y_grad_matrix(k,j) += ratio*y_grad[2][1];

      stiffness_matrix(i,k) += ratio*stiffness[0][2];
      stiffness_matrix(j,k) += ratio*stiffness[1][2];
      stiffness_matrix(k,i) += ratio*stiffness[2][0];
      stiffness_matrix(k,j) += ratio*stiffness[2][1];

      for (unsigned int l=0; l<dof_per_cell; ++l)
      {
        mass_matrix(k,l) += square_ratio*mass[2][2];
        x_grad_matrix(k,l) += square_ratio*x_grad[2][2];
        y_grad_matrix(k,l) += square_ratio*y_grad[2][2];
        stiffness_matrix(k,l) += square_ratio*stiffness[2][2];
      }
    }
  }
}
