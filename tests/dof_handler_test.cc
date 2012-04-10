#include <cassert>
#include <string>
#include "CELL.hh"
#include "DOF_HANDLER.hh"
#include "EDGE.hh"
#include "GLC.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
#include "TRIANGULATION.hh"

using namespace std;

typedef vector<unsigned int> ui_vector;

int main(int argc,char** argv)
{
  string geometry_inp("/home/bruno/Documents/Transport/janus/tests/geometry_dof.inp");
  string parameters_inp("/home/bruno/Documents/Transport/janus/tests/parameters_dof.inp");

  TRIANGULATION triangulation(&geometry_inp);
  PARAMETERS parameters(&parameters_inp);
  vector<QUADRATURE*> quad(2,NULL);
   
  triangulation.Read_geometry();
  triangulation.Build_edges();

  parameters.Read_parameters(triangulation.Get_n_sources(),
      triangulation.Get_n_materials());
  
  quad[0] = new GLC(parameters.Get_sn_order(),parameters.Get_L_max(),
      parameters.Get_galerkin());
  quad[0]->Build_quadrature();
  
  quad[1] = new GLC(parameters.Get_sn_order()/2.,parameters.Get_L_max()/2.,
      parameters.Get_galerkin());
  quad[1]->Build_quadrature();

  DOF_HANDLER dof_handler(&triangulation,parameters);
  dof_handler.Compute_sweep_ordering(quad);
 
  // Check the most normal direction for the bottom
  for (unsigned int i=0; i<12; ++i)
  {
    if (i==1 || i==10)    
      assert(dof_handler.Is_most_normal_bottom(0,i)==true);
    else
      assert(dof_handler.Is_most_normal_bottom(0,i)==false);
  }
  for (unsigned int i=0; i<4; ++i)
  {
    if (i==0 || i==3)    
      assert(dof_handler.Is_most_normal_bottom(1,i)==true);
    else
      assert(dof_handler.Is_most_normal_bottom(1,i)==false);
  }
  // Check the most normal direction for the right
  for (unsigned int i=0; i<12; ++i)
  {
    if (i==6 || i==9)    
      assert(dof_handler.Is_most_normal_right(0,i)==true);
    else
      assert(dof_handler.Is_most_normal_right(0,i)==false);
  }
  for (unsigned int i=0; i<4; ++i)
  {
    if (i==2 || i==3)    
      assert(dof_handler.Is_most_normal_right(1,i)==true);
    else
      assert(dof_handler.Is_most_normal_right(1,i)==false);
  }
  // Check the most normal direction for the top
  for (unsigned int i=0; i<12; ++i)
  {
    if (i==7 || i==4)    
      assert(dof_handler.Is_most_normal_top(0,i)==true);
    else
      assert(dof_handler.Is_most_normal_top(0,i)==false);
  }
  for (unsigned int i=0; i<4; ++i)
  {
    if (i==1 || i==2)    
      assert(dof_handler.Is_most_normal_top(1,i)==true);
    else
      assert(dof_handler.Is_most_normal_top(1,i)==false);
  }
  // Check the most normal direction for the right
  for (unsigned int i=0; i<12; ++i)
  {
    if (i==0 || i==3)    
      assert(dof_handler.Is_most_normal_left(0,i)==true);
    else
      assert(dof_handler.Is_most_normal_left(0,i)==false);
  }
  for (unsigned int i=0; i<4; ++i)
  {
    if (i==0 || i==1)    
      assert(dof_handler.Is_most_normal_left(1,i)==true);
    else
      assert(dof_handler.Is_most_normal_left(1,i)==false);
  }

  // Check the number of significant angular fluxes per direction
  assert(dof_handler.Get_n_sf_per_dir()==14);

  // Check the number of degrees of freedom
  assert(dof_handler.Get_n_dof()==20);

  // Check the number of cells
  assert(dof_handler.Get_n_cells()==5);

  // Chek the significan angular flux map
  assert(dof_handler.Get_saf_map_dof(0)==0);
  assert(dof_handler.Get_saf_map_dof(2)==6);
  assert(dof_handler.Get_saf_map_dof(4)==12);

  // Check the sweep ordering
  ui_vector const* sweep_ordering(dof_handler.Get_sweep_order(0,0));
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==2);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==1);
  assert((*sweep_ordering)[4]==3);

  sweep_ordering = dof_handler.Get_sweep_order(0,1);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==1);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==2);
  assert((*sweep_ordering)[4]==3);

  sweep_ordering = dof_handler.Get_sweep_order(0,2);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==2);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==1);
  assert((*sweep_ordering)[4]==3);

  sweep_ordering = dof_handler.Get_sweep_order(1,0);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==2);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==1);
  assert((*sweep_ordering)[4]==3);

  sweep_ordering = dof_handler.Get_sweep_order(0,3);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==2);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==3);
  assert((*sweep_ordering)[4]==1);

  sweep_ordering = dof_handler.Get_sweep_order(0,4);
  assert((*sweep_ordering)[0]==2);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==0);
  assert((*sweep_ordering)[4]==1);

  sweep_ordering = dof_handler.Get_sweep_order(0,5);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==2);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==3);
  assert((*sweep_ordering)[4]==1);

  sweep_ordering = dof_handler.Get_sweep_order(1,1);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==2);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==3);
  assert((*sweep_ordering)[4]==1);

  sweep_ordering = dof_handler.Get_sweep_order(0,6);
  assert((*sweep_ordering)[0]==1);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==2);
  assert((*sweep_ordering)[4]==0);

  sweep_ordering = dof_handler.Get_sweep_order(0,7);
  assert((*sweep_ordering)[0]==2);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==1);
  assert((*sweep_ordering)[4]==0);

  sweep_ordering = dof_handler.Get_sweep_order(0,8);
  assert((*sweep_ordering)[0]==1);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==2);
  assert((*sweep_ordering)[4]==0);

  sweep_ordering = dof_handler.Get_sweep_order(1,2);
  assert((*sweep_ordering)[0]==1);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==2);
  assert((*sweep_ordering)[4]==0);

  sweep_ordering = dof_handler.Get_sweep_order(0,9);
  assert((*sweep_ordering)[0]==1);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==0);
  assert((*sweep_ordering)[4]==2);

  sweep_ordering = dof_handler.Get_sweep_order(0,10);
  assert((*sweep_ordering)[0]==0);
  assert((*sweep_ordering)[1]==1);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==3);
  assert((*sweep_ordering)[4]==2);

  sweep_ordering = dof_handler.Get_sweep_order(0,11);
  assert((*sweep_ordering)[0]==1);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==0);
  assert((*sweep_ordering)[4]==2);

  sweep_ordering = dof_handler.Get_sweep_order(1,3);
  assert((*sweep_ordering)[0]==1);
  assert((*sweep_ordering)[1]==3);
  assert((*sweep_ordering)[2]==4);
  assert((*sweep_ordering)[3]==0);
  assert((*sweep_ordering)[4]==2);

  // Check the cell iterator
  vector<CELL*>::iterator cell(dof_handler.Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler.Get_mesh_end());
  for (unsigned int id=0; cell<cell_end; ++cell,++id)
   assert((*cell)->Get_id()==id); 

  // Check the edge iterator
  vector<EDGE>::iterator edge(dof_handler.Get_edges_begin());
  vector<EDGE>::iterator edge_end(dof_handler.Get_edges_end());
  for (unsigned int gid=0; edge<edge_end; ++edge,++gid)
    assert(edge->Get_gid()==gid);

  return 0;
}
