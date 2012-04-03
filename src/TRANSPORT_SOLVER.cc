#include "TRANSPORT_SOLVER.hh"

TRANSPORT_SOLVER::TRANSPORT_SOLVER(string* g_inputfile,string* p_inputfile,
    string* outputfile,Epetra_MpiComm* mpi_comm) :
  init_timer(NULL),
  calc_timer(NULL),
  parameters(p_inputfile),
  outputfile(outputfile),
  comm(mpi_comm),
  flux_moments(NULL),
  flux_moments_map(NULL),
  triangulation(g_inputfile),
  dof_handler(NULL)
{
  // Initialize the timers
  init_timer = new Teuchos::Time("Initialization time");
  calc_timer = new Teuchos::Time("Calculation time");
  // Start the init_timer
  init_timer->start();
  
  // Build the triangulation
  cout<<"Building the triangulation"<<endl;
  triangulation.Read_geometry();
  triangulation.Build_edges();

  // Read the parameters
  cout<<"Reading the parameters"<<endl;
  parameters.Read_parameters(triangulation.Get_n_sources(),
      triangulation.Get_n_materials());

  // Build the quadratures
  cout<<"Building the quadratures"<<endl;
  unsigned int n_lvl(parameters.Get_n_levels());
  unsigned int tmp_sn(parameters.Get_sn_order());
  unsigned int tmp_L_max(parameters.Get_L_max());
  if (parameters.Get_mip()==true || parameters.Get_multigrid()==true)
    --n_lvl;
  quad.resize(n_lvl,NULL);
  for (unsigned int lvl=0; lvl<n_lvl; ++lvl)
  {
    if (parameters.Get_quad_type()==ls)
      quad[lvl] = new LS(tmp_sn,tmp_L_max,parameters.Get_galerkin());
    else
      quad[lvl] = new GLC(tmp_sn,tmp_L_max,parameters.Get_galerkin());
    quad[lvl]->Build_quadrature();

    tmp_sn = ceil(tmp_sn/2.);
    if (tmp_L_max>tmp_sn)
      tmp_L_max = tmp_sn;
  }                               

  // Build the dof handler
  cout<<"Building the dof handler"<<endl;
  dof_handler = new DOF_HANDLER(&triangulation,parameters);
  cout<<"Compute sweep ordering"<<endl;
  dof_handler->Compute_sweep_ordering(quad);


  // Instantiate the flux moments map and vector
  cout<<"Instantiate the flux moments map and vector"<<endl;
  flux_moments_size = dof_handler->Get_n_dof()*quad[0]->Get_n_mom()+
    dof_handler->Get_n_sf_per_dir()*quad[0]->Get_n_dir();
  flux_moments_map = new Epetra_Map(flux_moments_size,0,*comm);
  flux_moments = new Epetra_MultiVector(*flux_moments_map,1);
  cout<<"Initialization done"<<endl;
  init_timer->stop();
}

TRANSPORT_SOLVER::~TRANSPORT_SOLVER()
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

void TRANSPORT_SOLVER::Solve()
{
  calc_timer->start();
  double building_mip_time(0.);
  double solve_mip_time(0.);
  // Compute rhs
  if (parameters.Get_multigrid()==true)
  {
  }
  else
  {
    if (parameters.Get_solver_type()!=si)
    {
      TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,quad[0],
          comm,flux_moments_map);

      // Compute right-hand side of GMRES (uncollided flux moments (S*inv(L)*q))
      Epetra_MultiVector rhs(*flux_moments);
      transport_operator.Sweep(rhs,true);

      Epetra_LinearProblem problem(&transport_operator,flux_moments,&rhs);

      AztecOO solver(problem);

      // By default use ILUT as preconditioner -> by default the Krylov
      // solvers cannot be used matrix-free
      solver.SetAztecOption(AZ_precond,AZ_none);

      if (parameters.Get_solver_type()==bicgstab)
        solver.SetAztecOption(AZ_solver,AZ_bicgstab);
      else 
      {
        // Restart parameter for GMRES
        int krylov_space(30);
        // Set a new Krylov subspace size
        solver.SetAztecOption(AZ_kspace,krylov_space);
        if (parameters.Get_solver_type()==gmres_condnum)
          solver.SetAztecOption(AZ_solver,AZ_gmres_condnum);
        else
          solver.SetAztecOption(AZ_solver,AZ_gmres);
      }
      // Convergence criterion ||r||_2/||b||_2
      solver.SetAztecOption(AZ_conv,AZ_rhs);
      
      // Solve the transport equation
      solver.Iterate(parameters.Get_max_it(),parameters.Get_tolerance());

      if (parameters.Get_mip()==true)
      {
        MIP precond(dof_handler,&parameters,quad[0],comm);
        precond.Solve(*flux_moments);
        building_mip_time = precond.Get_building_mip_time();
        solve_mip_time = precond.Get_solve_mip_time();
      }
    }
    else
    {
      const unsigned int max_it(parameters.Get_max_it());
      double convergence(1.);

      TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,quad[0],
          comm,flux_moments_map);

      Epetra_MultiVector flux_moments_old(*flux_moments);

      for (unsigned int i=0; i<max_it; ++i)
      {
        double num(1.);
        double denom(1.);

        transport_operator.Compute_scattering_source(*flux_moments);
        transport_operator.Sweep(*flux_moments,true);

        Epetra_MultiVector diff_flux(*flux_moments);
        Epetra_BLAS blas;
        blas.AXPY(flux_moments_size,-1.,flux_moments_old.Values(),diff_flux.Values());

        // To change: use the norm of the residual for the convergence i.e.
        // the same criterion as the krylov vector
        assert(diff_flux.Norm2(&num)==0);
        assert(flux_moments->Norm2(&denom)==0);
        convergence = num/denom;
        cout<<"Convergence at iteration "<<i<<": "<<convergence<<endl;
        if (convergence<parameters.Get_tolerance())
          break;

        flux_moments_old = *flux_moments;
      }
    }
  }
  calc_timer->stop();
  cout<<"Initialization time: "<<init_timer->totalElapsedTime()<<endl;
  cout<<"Calculation time: "<<calc_timer->totalElapsedTime()<<endl;
  if (parameters.Get_mip()==true)
  {
    cout<<"Building MIP time: "<<building_mip_time<<endl;
    cout<<"Solving MIP time: "<<solve_mip_time<<endl;
  }
  cout<<"Total elapsed time: "<<init_timer->totalElapsedTime()+
    calc_timer->totalElapsedTime()<<endl;
}

void TRANSPORT_SOLVER::Write_in_file()
{
  const unsigned int n_cells(dof_handler->Get_n_cells());
  const unsigned int n_dof(dof_handler->Get_n_dof());
  const unsigned int n_mom(quad[0]->Get_n_mom());
  d_vector offset(n_cells+1,0.);
  vector<d_vector> points;
  vector<CELL*>::iterator cell_it(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());

  for (; cell_it<cell_end; ++cell_it)
  {
    const unsigned int cell_id((*cell_it)->Get_id());
    FINITE_ELEMENT const* const fe((*cell_it)->Get_fe());
    vector<d_vector> tmp_points((*cell_it)->Reorder_vertices());
    const unsigned int tmp_pts_size(tmp_points.size());
    for (unsigned int i=0; i<tmp_pts_size; ++i)
      points.push_back(tmp_points[i]);
    if (cell_id<n_cells-1)
      offset[cell_id+1] = offset[cell_id]+fe->Get_dof_per_cell();
  }
  offset[n_cells] = n_dof;

  // Write the mesh and the flux in a text file
  string filename("silo_");
  filename += *outputfile;
  filename += ".txt";
  ofstream file(filename.c_str(),ios::out | ios::trunc);

  file<<n_cells<<" ";
  file<<n_dof<<" ";
  file<<quad[0]->Get_n_mom()<<"\n";

  for (unsigned int i=0; i<=n_cells; ++i)
    file<<offset[i]<<" ";
  file<<"\n";

  const unsigned int n_points(points.size());
  for (unsigned int i=0; i<n_points; ++i)
    file<<points[i][0]<<" "<<points[i][1]<<"\n";

  const unsigned int n_values(n_dof*n_mom);
  for (unsigned int i=0; i<n_values; ++i)
    file<<(*flux_moments)[0][i]<<"\n";

  // To have the scalar flux phi_00 needs to be multiply by sqrt(4*PI)
  for (unsigned int i=0; i<n_dof; ++i)
    file<<(*flux_moments)[0][i]*2.*M_SQRTPI<<"\n";

  file.close();
}
