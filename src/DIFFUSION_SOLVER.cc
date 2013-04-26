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

#include "DIFFUSION_SOLVER.hh"

DIFFUSION_SOLVER::DIFFUSION_SOLVER(string* g_inputfile,string* p_inputfile,
    string* xs_inputfile,string* outputfile,Epetra_MpiComm* mpi_comm) :
  init_timer(NULL),
  calc_timer(NULL),
  outputfile(outputfile),
  comm(mpi_comm),
  group_flux(NULL),
  scalar_flux_map(NULL),
  cross_sections(xs_inputfile),
  parameters(p_inputfile),
  triangulation(g_inputfile),
  dof_handler(NULL)
{
  // Initialize the timers
  init_timer = new Teuchos::Time("Initialization time");
  calc_timer = new Teuchos::Time("Calculation time");
  // Start the init_timer
  init_timer->start();
  
  // Build the triangulation
  cout<<"Start building the triangulation"<<endl;
  triangulation.Read_geometry();
  triangulation.Build_edges();
  cout<<"Done building the triangulation"<<endl;

  // Read the parameters
  cout<<"Start initializing the parameters"<<endl;
  parameters.Read_diffusion_parameters(triangulation.Get_n_sources());
  cout<<"Done initializing the parameters"<<endl;

  // Read the cross sections
  cout<<"Start reading the cross sections"<<endl;
  PERMUTATION_TYPE permutation_type(none);
  cross_sections.Read_regular_cross_sections(triangulation.Get_n_materials(),
      permutation_type,false);
  cross_sections.Apply_ang_lvls_and_tc(parameters.Get_multigrid(),
      parameters.Get_transport_correction(),parameters.Get_optimal_tc(),
      triangulation.Get_n_materials(),parameters.Get_sn_order());
  cout<<"Done reading the cross sections"<<endl;

  // Build the dof handler
  cout<<"Start building the dof handler"<<endl;
  dof_handler = new DOF_HANDLER(&triangulation,parameters,cross_sections);
  cout<<"Done building the dof handler"<<endl;


  // Instantiate the scalar flux  map and group_flux
  scalar_flux_size = dof_handler->Get_n_dof();
  scalar_flux_map = new Epetra_Map(scalar_flux_size,0,*comm);
  group_flux = new Epetra_MultiVector(*scalar_flux_map,cross_sections.Get_n_groups());
  init_timer->stop();
}

DIFFUSION_SOLVER::~DIFFUSION_SOLVER()
{
  if (group_flux!=NULL)
  {
    delete group_flux;
    group_flux = NULL;
  }
  
  if (scalar_flux_map!=NULL)
  {
    delete scalar_flux_map;
    scalar_flux_map = NULL;
  }

  if (dof_handler!=NULL)
  {
    delete dof_handler;
    dof_handler = NULL;
  }

  if (calc_timer!=NULL)
  {
    delete calc_timer;
    calc_timer = NULL;
  }
  
  if (init_timer!=NULL)
  {
    delete calc_timer;
    calc_timer = NULL;
  }
}

void DIFFUSION_SOLVER::Solve()
{
  calc_timer->start();
  unsigned int group_iter(0);
  const unsigned int max_group_it(parameters.Get_max_group_it());
  const unsigned int n_groups(cross_sections.Get_n_groups());
  const unsigned int n_refinements(parameters.Get_n_refinements());
  const double group_tol(parameters.Get_group_tolerance());
  for (unsigned int r=0; r<=n_refinements; ++r)
  {
    const unsigned int n_dof(dof_handler->Get_n_dof());
    double group_conv(10.*group_tol);
    Epetra_MultiVector old_group_flux(*group_flux);
    Epetra_MultiVector scalar_flux(*scalar_flux_map,1);
    Epetra_MultiVector initial_guess(*scalar_flux_map,1);

    MIP mip(comm,&parameters,dof_handler);
    // Loop over the groups
    while (group_conv>group_tol)
    {
      // Loop over the supergroups.
      for (unsigned int g=0; g<n_groups; ++g)
      {
        // Set the current group
        mip.Set_group(g);
        //copy((*group_flux)[g],(*group_flux)[g]+n_dof,initial_guess[0]);
        //mip.Solve_diffusion(n_groups,scalar_flux,*group_flux,&initial_guess);
        mip.Solve_diffusion(n_groups,scalar_flux,*group_flux);

        copy((*group_flux)[g],(*group_flux)[g]+n_dof,old_group_flux[g]);
        copy(scalar_flux[0],scalar_flux[0]+n_dof,(*group_flux)[g]);
      }
      // Compute the convergence over all the groups
      group_conv = Compute_convergence(*group_flux,old_group_flux,n_groups);
      cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
        group_conv<<endl;
      ++group_iter;
      if (group_iter>=max_group_it)
        break;
    }
    // Refine the mesh
    if (r<n_refinements)
      Refine_mesh(r);
  }
  calc_timer->stop();
  cout<<"Initialization time: "<<init_timer->totalElapsedTime()<<endl;
  cout<<"Calculation time: "<<calc_timer->totalElapsedTime()<<endl;
  cout<<"Total elapsed time: "<<init_timer->totalElapsedTime()+
    calc_timer->totalElapsedTime()<<endl;
}

void DIFFUSION_SOLVER::Refine_mesh(unsigned int r)
{
  ui_set cells_to_refine;
  ui_set adjacent_cells;
  vector<ui_vector> projection;
  map<unsigned int,vector<vector<d_vector> > > edge_to_refine;

  if (parameters.Get_verbose()>0)
  {
    stringstream r_string;
    r_string<<r;
    string filename("silo_flux_r"+r_string.str()+"_before.txt");    
    Write_in_file(&filename);
  }

  // Compute the error estimate and flag the cells
  ERROR_ESTIMATOR::Compute_refinement(cross_sections.Get_n_groups(),
      dof_handler,&parameters,group_flux,cells_to_refine,adjacent_cells,
      edge_to_refine);
  // Build the refined grid
  triangulation.Refine_mesh(cells_to_refine,adjacent_cells,edge_to_refine,
      projection);
  // Project the solution 
  Project_solution(projection);
  // Free the dof_handler
  delete dof_handler;
  // Build the new edges
  triangulation.Build_edges();
  // Create the new dof_handler
  dof_handler = new DOF_HANDLER(&triangulation,parameters,cross_sections);

  if (parameters.Get_verbose()>0)
  {
    stringstream r_string;
    r_string<<r;
    string filename("silo_flux_r"+r_string.str()+"_after.txt");    
    Write_in_file(&filename);
  }
}

void DIFFUSION_SOLVER::Project_solution(vector<ui_vector> const &projection)
{
  scalar_flux_size = projection.size();
  Epetra_Map* new_scalar_flux_map = new Epetra_Map(scalar_flux_size,0,*comm);
  Epetra_MultiVector* new_group_flux = new Epetra_MultiVector(*new_scalar_flux_map,
      cross_sections.Get_n_groups());

 // // Project group_flux on the new grid
 // for (unsigned int g=0; g<cross_sections.Get_n_groups(); ++g)
 // {
 //   for (unsigned int i=0; i<scalar_flux_size; ++i)
 //   {
 //     unsigned int n_old_vertices(projection[i].size());
 //     for (unsigned int j=0; j<n_old_vertices; ++j)
 //       (*new_group_flux)[g][i] += (*group_flux)[g][projection[i][j]];
 //     (*new_group_flux)[g][i] /= (double)(n_old_vertices);
 //   }
 // }

  // Copy new_scalar_flux_map and new_group_flux to scalar_flux_map and
  // group_flux
  delete group_flux;
  delete scalar_flux_map;
  scalar_flux_map = new_scalar_flux_map;
  group_flux = new_group_flux;
  new_scalar_flux_map = NULL;
  new_group_flux = NULL;
}

double DIFFUSION_SOLVER::Compute_convergence(Epetra_MultiVector const &flux, 
    Epetra_MultiVector const &old_flux, const unsigned int n) const
{
  double conv(0.);
  double* l2_norm_num = new double[n];
  double* l2_norm_denom = new double[n];
  Epetra_MultiVector tmp(flux);
  tmp.Update(-1.,old_flux,0.);
  // Compute the L2 norm of each vector of the MultiVector
  tmp.Norm2(l2_norm_num);
  flux.Norm2(l2_norm_denom);
  double norm_num(0.);
  double norm_denom(0.);
  for (unsigned int j=0; j<n; ++j)
  {
    norm_num += pow(l2_norm_num[j],2);
    norm_denom += pow(l2_norm_denom[j],2);
  }
  conv = sqrt(norm_num/norm_denom);
  delete [] l2_norm_num;
  delete [] l2_norm_denom;      

  return conv;
}

void DIFFUSION_SOLVER::Write_in_file(string* filename_) const
{
  const unsigned int n_cells(dof_handler->Get_n_cells());
  const unsigned int n_dof(dof_handler->Get_n_dof());
  const unsigned int n_groups(cross_sections.Get_n_groups());
  ui_vector offset(n_cells+1,0.);
  vector<d_vector> points;
  vector<d_vector> c_points;
  points.reserve(n_dof);
  c_points.reserve(n_dof+n_cells);
  vector<CELL*>::iterator cell_it(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
  Epetra_Map c_map(n_cells,0,*comm);
  Epetra_MultiVector c_flux(c_map,cross_sections.Get_n_groups());

  for (; cell_it<cell_end; ++cell_it)
  {
    const unsigned int cell_id((*cell_it)->Get_id());
    FINITE_ELEMENT const* const fe((*cell_it)->Get_fe());
    vector<d_vector> tmp_points((*cell_it)->Reorder_vertices());
    d_vector center(2,0.);
    const unsigned int tmp_pts_size(tmp_points.size());
    for (unsigned int i=0; i<tmp_pts_size; ++i)
    {
      points.push_back(tmp_points[i]);
      center[0] += tmp_points[i][0]/static_cast<const double> (tmp_pts_size);
      center[1] += tmp_points[i][1]/static_cast<const double> (tmp_pts_size);
    }
    c_points.push_back(center);

    if (cell_id<n_cells-1)
      offset[cell_id+1] = offset[cell_id]+fe->Get_dof_per_cell();
  }
  offset[n_cells] = n_dof;

  // Write the mesh and the flux in a text file
  string filename;
  if (filename_==NULL)
    filename = "silo_"+ *outputfile+".txt";
  else
    filename = *filename_;
  ofstream file(filename.c_str(),ios::out | ios::trunc);

  file<<n_cells<<" ";
  file<<n_dof<<" ";
  file<<n_groups<<" ";

  for (unsigned int i=0; i<=n_cells; ++i)
    file<<offset[i]<<" ";
  file<<"\n";

  const unsigned int n_points(points.size());
  for (unsigned int i=0; i<n_points; ++i)
    file<<points[i][0]<<" "<<points[i][1]<<"\n";

  // Loop over the group
  for (unsigned int g=0; g<n_groups; ++g)
    for (unsigned int i=0; i<n_dof; ++i)
      file<<(*group_flux)[g][i]<<"\n";

  const unsigned int n_cpoints(c_points.size());
  for (unsigned int i=0; i<n_cpoints; ++i)
    file<<c_points[i][0]<<" "<<c_points[i][1]<<"\n";

  for (unsigned int g=0; g<n_groups; ++g)
  {
    for (unsigned int i=0; i<n_cells; ++i)
    {
      double center_flux(0.0);
      double n_vertices(offset[i+1]-offset[i]);
      for (unsigned int j=offset[i]; j<offset[i+1]; ++j)
        center_flux += (*group_flux)[g][j];
      c_flux[g][i] = center_flux/n_vertices;
    }
  }
  for (unsigned int g=0; g<n_groups; ++g)
    for (unsigned int i=0; i<n_cpoints; ++i)
      file<<c_flux[g][i]<<"\n";

  file.close();
}
