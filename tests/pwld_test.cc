#include <cassert>
#include <cmath>
#include <vector>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "CELL.hh"
#include "PWLD.hh"

using namespace std;

typedef vector<double> d_vector;

int main(int argc,char** argv)
{
  unsigned int dof_per_cell(4);
  double m_r(1./96.);
  double g_r(1./12.);
  double x0,x1,x2,y0,y1,y2,x1_x0,y2_y0,y1_y0,x2_x0;
  double a_00,a_01,a_01,a_02,a_11,a_22,jacobian;
  d_vector cell_x(4,0.);
  d_vector cell_y(4,0.);
  vector<d_vector> mass_matrix(4,d_vector(4,0.));
  vector<d_vector> x_grad_matrix(4,d_vector(4,0.));
  vector<d_vector> y_grad_matrix(4,d_vector(4,0.));
  vector<d_vector> sub_stiffness_matrix(3,d_vector(3,0.));
  vector<d_vector> stiffness_matrix(4,d_vector(4,0.));
  cell_x[0] = 1.;
  cell_x[1] = 2.;
  cell_x[2] = 2.;
  cell_x[3] = 1.;
  cell_y[0] = 2.;
  cell_y[1] = 2.;
  cell_y[2] = 3.; 
  cell_y[3] = 3.;
  x0 = cell_x[0];
  x1 = cell_x[1];
  x2 = (cell_x[0]+cell_x[1]+cell_x[2]+cell_x[3])/4.;
  y0 = cell_y[0];
  y1 = cell_y[1];
  y2 = (cell_y[0]+cell_y[1]+cell_y[2]+cell_y[3])/4.;
  x1_x0 = x1-x0;
  x2_x0 = x2-x0;
  y2_y0 = y2-y0;
  y1_y0 = y1-y0;
  a_00(0.5*(pow(y2_y1,2)+pow(x2_x1,2)));
  a_01(-0.5*(y2_y0*y2_y1+x2_x0*x2_x1)); 
  a_02(0.5*(y1_y0*y2_y1+x1_x0*x2_x1));  
  a_11(0.5*(pow(y2_y0,2)+pow(x2_x0,2)));
  a_12(-0.5*(y2_y0*y1_y0+x2_x0*x1_x0)); 
  a_22(0.5*(pow(y1_y0,2)+pow(x1_x0,2)));
  jacobian = fabs(x1_x0*y2_y0-y1_y0*x2_x0);
  stiffness[0][0] = a_00/jacobian;
  stiffness[0][1] = a_01/jacobian;
  stiffness[0][2] = a_02/jacobian;
  stiffness[1][0] = a_01/jacobian;
  stiffness[1][1] = a_11/jacobian;
  stiffness[1][2] = a_12/jacobian;
  stiffness[2][0] = a_02/jacobian;
  stiffness[2][1] = a_12/jacobian;
  stiffness[2][2] = a_22/jacobian;
  PWLD fe(cell_x,cell_y);
  mass_matrix[0][0] = 11.*m_r;
  mass_matrix[0][1] = 5.*m_r;
  mass_matrix[0][2] = 3.*m_r;
  mass_matrix[0][3] = 5.*m_r;
  mass_matrix[1][0] = 5.*m_r;
  mass_matrix[1][1] = 11.*m_r;
  mass_matrix[1][2] = 5.*m_r;
  mass_matrix[1][3] = 3.*m_r;
  mass_matrix[2][0] = 3.*m_r;
  mass_matrix[2][1] = 5.*m_r;
  mass_matrix[2][2] = 11.*m_r;
  mass_matrix[2][3] = 5.*m_r;
  mass_matrix[3][0] = 5.*m_r;
  mass_matrix[3][1] = 3.*m_r;
  mass_matrix[3][2] = 5.*m_r;
  mass_matrix[3][3] = 11.*m_r;
  x_grad_matrix[0][0] = -2.*g_r;
  x_grad_matrix[0][1] = -2.*g_r;
  x_grad_matrix[0][2] = -1.*g_r;
  x_grad_matrix[0][3] = -1.*g_r;
  x_grad_matrix[1][0] = 2.*g_r;
  x_grad_matrix[1][1] = 2.*g_r;
  x_grad_matrix[1][2] = 1.*g_r;
  x_grad_matrix[1][3] = 1.*g_r;
  x_grad_matrix[2][0] = 1.*g_r;
  x_grad_matrix[2][1] = 1.*g_r;
  x_grad_matrix[2][2] = 2.*g_r;
  x_grad_matrix[2][3] = 2.*g_r;
  x_grad_matrix[3][0] = -1.*g_r;
  x_grad_matrix[3][1] = -1.*g_r;
  x_grad_matrix[3][2] = -2.*g_r;
  x_grad_matrix[3][3] = -2.*g_r;
  y_grad_matrix[0][0] = -2.*g_r;
  y_grad_matrix[0][1] = -1.*g_r;
  y_grad_matrix[0][2] = -1.*g_r;
  y_grad_matrix[0][3] = -2.*g_r;
  y_grad_matrix[1][0] = -1.*g_r;
  y_grad_matrix[1][1] = -2.*g_r;
  y_grad_matrix[1][2] = -2.*g_r;
  y_grad_matrix[1][3] = -1.*g_r;
  y_grad_matrix[2][0] = 1.*g_r;
  y_grad_matrix[2][1] = 2.*g_r;
  y_grad_matrix[2][2] = 2.*g_r;
  y_grad_matrix[2][3] = 1.*g_r;
  y_grad_matrix[3][0] = 2.*g_r;
  y_grad_matrix[3][1] = 1.*g_r;
  y_grad_matrix[3][2] = 1.*g_r;
  y_grad_matrix[3][3] = 2.*g_r;

  // First sub-triangle
  stiffness_matrix[0][0] += a_00 + 1./4.*a_02 + 1./4.*a_02 + 1./16.*a_22;
  stiffness_matrix[0][1] += a_01 + 1./4.*a_02 + 1./4.*a_12 + 1./16.*a_22;

  cout<<"before fe_1d"<<endl;
  fe.Build_fe_1d();
  //fe.Build_upwind_matrices(cell,mesh);
  cout<<"before fe_2d"<<endl;
  fe.Build_fe_2d();
  cout<<"after fe_2d"<<endl;

  // Check the number of degrees of freedom per cell
  assert(dof_per_cell==fe.Get_dof_per_cell());

  // Check the mass_matrix
  Teuchos::SerialDenseMatrix<int,double> const* mass(fe.Get_mass_matrix());
  for (unsigned int i=0; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      assert(fabs(mass_matrix[i][j]-(*mass)(i,j))<1e-12);
  // Check the x_grad_matrix
  Teuchos::SerialDenseMatrix<int,double> const* x_grad(fe.Get_x_grad_matrix());
  for (unsigned int i=0; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      assert(fabs(x_grad_matrix[i][j]-(*x_grad)(i,j))<1e-12);
  // Check the y_grad_matrix
  Teuchos::SerialDenseMatrix<int,double> const* y_grad(fe.Get_y_grad_matrix());
  for (unsigned int i=0; i<4; ++i)
    for (unsigned int j=0; j<4; ++j)
      assert(fabs(y_grad_matrix[i][j]-(*y_grad)(i,j))<1e-12);
  // Check the stiffness_matrix
  // TO DO

  return 0;
}
