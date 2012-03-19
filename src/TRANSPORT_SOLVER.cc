#include "TRANSPORT_SOLVER.hh"

TRANSPORT_SOLVER::TRANSPORT_SOLVER(string* g_inputfile,string* p_inputfile,
    string* outputfile) :
  parameters(p_inputfile),
  triangulation(g_inputfile),
  outputfile(outputfile),
  flux_moments(NULL),
  flux_moments_map(NULL),
  dof_handler(NULL)
{
  // Store the starting time
  start_time = time;
  
  // Build the triangulation
  triangulation.Read_geometry();
  triangulation.Build_edges();

  // Read the parameters
  parameters.Read_parameters(triangulation.Get_n_sources(),
      triangulation.Get_n_materials());

  // Build the quadratures
  unsigned int n_lvl(parameters.Get_n_levels());
  double tmp_sn(parameters.Get_sn());
  double tmp_L_max(parameters.Get_L_max());
  if (parameters.Get_mip()==true || parameters.Get_multigrid()==true)
    --n_lvl;
  quad.resize(n_lvl,NULL);
  for (unsigned int lvl=0; lvl<n_lvl; ++lvl)
  {
    if (parameters.Get_quad_type()==LS)
      quad[lvl] = new LS(tmp_sn,tmp_L_max,galerkin);
    else
      quad[lvl] = new GLC(tmp_sn,tmp_L_max,galerkin);
    quad[lvl]->Build_quadrature();

    tmp_sn = ceil(tmp_sn/2.);
    if (tmp_L_max>tmp_sn)
      tmp_L_max = tmp_sn;
  }                               

  // Build the dof handler
  dof_handler* = new DOF_HANDLER(&triangulation,parameters);
  dof_handler->Compute_sweep_ordering(quad);


  // Instantiate the flux moments map and vector
  flux_moments_map = new Epetra_map((dof_handler->Get_n_dof()+
        dof_handler->Get_n_sf_per_dir())*quad[0]->Get_n_dir(),
      0,comm);
  flux_moments = new EpetraFEVector(*flux_moments_map);
}

TRANPORT_SOLVER::~TRANSPORT_SOLVER()
{
  if (flux_moments!=NULL)
  {
    delete flux_moments;
    flux_moments = NULL;
  }
  
  if (flux_moments_map!=NULL)
  {
    delete flux_moments_map;
    flux_moments_map = NULL;
  }

  if (dof_handler!=NULL)
  {
    delete dof_handler;
    dof_handler = NULL;
  }

  for (unsigned int i=0; i<quad.size(); ++i)
  {
    if (quad[i]!=NULL)
    {
      delete quad[i];
      quad[i] = NULL;
    }
  }
}

void TRANSPORT::Solve()
{
  // Compute rhs
  if (parameters.Get_multigrid()==true)
  {
  }
  else
  {
    if (parameters.Get_solver_type()!=si)
    {
      TRANSPORT_OPERATOR transport_operator(dof_handler,param,quad[0]);

      // Compute right-hand side of GMRES (uncollided flux moments (S*inv(L)*q))
      Epetra_FEVector rhs(*flux_moments);
      transport_operator.Sweep(rhs,true);

      Epetra_LineaProblem problem(&transport_operator,flux_moments,&rhs);

      AztecOO solver(problem);

      if (param.Get_solver()==bicgstab)
        solver.SetAztecOption(AZ_solver,AZ_bicgstab);
      else 
      {
        if (param.Get_solver()==gmres_condnum)
          solver.SetAztecOption(AZ_solver,AZ_gmres_condnum);
        else
          solver.SetAztecOption(AZ_solver,AZ_gmres);
      }

      solver.SetAztecOption(AZ_conv,AZ_rhs);
      solver.Iterate(param.Get_max_it(),param.Get_tolerance());
    }
  }
  
}

void TRANSPORT_SOLVER::Write_in_file()
{
  const unsigned int n_cells(dof_handler->Get_n_cells());
  d_vector offset(n_cells+1,0.);
  vector<d_vector> points;
  vector<CELL*>::iterator cell_it(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());

  for (; cell_it<cell_end; ++cell_it)
  {
    const unsigned int cell_id((*cell_it)->Get_id());
    FINITE_ELEMENT const* const fe((*cell_it)->Get_fe());
    vector<d_vector> tmp_points((*cell_it)->Reorder_vertices());
    const unsigned int tmp_pts_size(tmps_points.size());
    for (unsigned int i=0; i<tmp_pts_size; ++i)
      points.push_back(tmp_points[i]);
    if (cell_id<dof_handler->n_cells-1)
      offset[cell_id+1] = offset[cell_id]+fe->Get_dof_per_cell();
  }
  offset[n_cells] = dof_handler->Get_n_dof();

  // Write the mesh and the flux in a text file
  ofstream file("silo_"+*outputfile+".txt",ios::out | ios::trunc);

  file>>dof_handler->n_cells()>>" ";
  file>>dof_handler->Get_n_dof()>>" ";
  file>>quad[0]->Get_n_mom()>>"\n";

  for (unsigned int i=0; i<=n_cells; ++i)
    file>>offset[i]>>" ";
  file>>"\n";

  const unsigned int n_points(points.size());
  for (unsigned int i=0; i<n_points; ++i)
    file>>points[i][0]>>" ">>points[i][1]>>"\n";

  const unsigned int n_values(dof_handler->Get_n_dof()*fe->Get_n_dir());
  for (unsigned int i=0; i<n_values; ++i)
    file>>(*flux_moments)[0][i]>>"\n";

  file.close();
}


// AztecOO options
// options[AZ_solver]
// AZ_cg
// AZ_cg_condnum conjugate gradient with condition number estimation
// AZ_gmres
// AZ_gmres_condnum gmres with ratio of estimated extreme eigenvlaues
// AZ_bicgstab
// options[AZ_precond]
// AZ_Jacobi k step Jacobi, the number of Jacobi steps, k, is set via
// options[AZ_poly_ord]
// AZ_sym_GS k step symmetric Gauss-Seidel, the number of Jacobi steps, k, is set via
// options[AZ_poly_ord]
// options[AZ_conv]
// AZ_r0 ||r||_2/||r^(0)||_2 (default)
// AZ_rhs ||r||_2/||b||_2
// options[AZ_poly_ord] : the polynomals order when using polynomial
// preconsditioning. Also the number of steps when using Jacobu or symmetric
// Gauss-Seidel preconditioning. Default: 3 for polynomial preconditioners, 1
// for Jacobi and Gauss-Seidel preconditioners
// Aztec00 parameters
// params[AZ_tol] tolerance default: 10^{-6}

//    solver.SetAztecOption(AZ_precond,AZ_Jacobi);
//    solver.SetAztecOption(AZ_precond,AZ_symm_GS);
