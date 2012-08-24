#include "TRANSPORT_SOLVER.hh"

TRANSPORT_SOLVER::TRANSPORT_SOLVER(string* g_inputfile,string* p_inputfile,
    string& xs_inputfile,string* outputfile,Epetra_MpiComm* mpi_comm) :
  init_timer(NULL),
  calc_timer(NULL),
  outputfile(outputfile),
  comm(mpi_comm),
  flux_moments(NULL),
  flux_moments_map(NULL),
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
  parameters.Read_parameters(triangulation.Get_n_sources(),
      triangulation.Get_n_materials());
  cout<<"Done initializing the parameters"<<endl;

  // Read the cross sections
  cout<<"Start reading the cross sections"<<endl;
  switch (parameters.Get_xs_type())
  {
    case fp :
      {
        cross_sections.Build_fokker_planck_xs(triangulation.Get_n_materials());
        break;
      }
    case regular :
      {
        cross_sections.Read_regular_cross_sections(triangulation.Get_n_materials(),
            false);
        break;
      }
    case regular_exs :
      {
        cross_sections.Read_regular_cross_sections(triangulation.Get_n_materials(),
            true);
        break;
      }
    case cepxs :
      {         
        cross_sections.Read_cepxs_cross_sections(triangulation.Get_n_materials());
        break;
      }
    default :
      {
        Check(false,"Unknown type of cross section file.");
      }
  }
  cross_sections.Apply_ang_lvls_and_tc(parameters.Get_multigrid(),
      parameters.Get_transport_correction(),parameters.Get_optimal_tc());
  cout<<"Done reading the cross sections"<<endl;

  // Build the quadratures
  cout<<"Start building the quadratures"<<endl;
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
    quad[lvl]->Build_quadrature(parameters.Get_weight_sum());

    tmp_sn = ceil(tmp_sn/2.);
    if (tmp_L_max>tmp_sn)
      tmp_L_max = tmp_sn;
  }                               
  cout<<"Done building the quadratures"<<endl;

  // Build the dof handler
  cout<<"Start building the dof handler"<<endl;
  dof_handler = new DOF_HANDLER(&triangulation,parameters,cross_sections);
  dof_handler->Compute_sweep_ordering(quad);
  cout<<"Done building the dof handler"<<endl;


  // Instantiate the flux moments map and vector
  flux_moments_size = dof_handler->Get_n_dof()*quad[0]->Get_n_mom()+
    dof_handler->Get_n_sf_per_dir()*quad[0]->Get_n_dir();
  flux_moments_map = new Epetra_Map(flux_moments_size,0,*comm);
  flux_moments = new Epetra_MultiVector(*flux_moments_map,1);
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
  unsigned int group_iter(0);
  const unsigned int max_group_it(Get_max_group_it());
  const unsigned int max_supergroup_it(Get_max_supergroup_it());
  const unsigned int supergroup(parameters.Get_n_grps_in_supergrp());
  double building_mip_time(0.);
  double init_prec_mip_time(0.);
  double solve_mip_time(0.);
  const double inner_tol(parameters.Get_inner_tolerance());
  const double group_tol(parameters.Get_group_tolerance());
  // Compute rhs
  if (parameters.Get_multigrid()==true)
  {
    const unsigned int lvl(0);
    const unsigned int max_lvl(parameters.Get_n_levels()-1);
    double group_conv(10.*parameters.Get_group_tol());
    Epetra_MultiVector group_flux(flux_moments_map,
        parameters.Get_n_groups());
    Epetra_MultiVector old_group_flux(group_flux);

    TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,&quad,comm,
        flux_moments_map,lvl,max_lvl);
    while (group_conv>parameters.Get_group_tolerance())
    {
      // Loop over the supergroups.
      for (unsigned int i=0; i<parameters.Get_n_supergroups(); ++i)
      {
        unsigned int supergroup_iter(0);
        double supergroup_conv(10.*group_tol);
        Epetra_MultiVector supergroup_flux(flux_moments_map,
            parameters.Get_n_groups());
        Epetra_MultiVector old_supergroup_flux(group_flux_moments);
        while (supergroup_conv>group_tol)
        {   
          // Loop over the groups in a supergroup.
          for (unsigned int j=0; j<supergroup; ++j)
          {
            unsigned int g(i*supergroup+j);
            // Set the current group
            transport_operator.Set_group(g);
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
            solver.Iterate(parameters.Get_max_it(),inner_tol);

            // Store the new flux_moments
            old_supergroup_flux[j] = supergroup_flux[j];
            supergroup_flux[j] = flux_moments;
          }
          // Compute the convergence of the supergroup
          supergroup_conv = Compute_convergence(supergroup_flux,
              old_supergroup_flux,supergroup);
          if (parameters.Get_verbose()>0)
          cout<<"Convergence over the supergroup "<<i<<" at iteration "<<
            supergroup_iter<<": "<<supergroup_conv<<endl;
          if (supergroup_iter>=max_supergroup_it)
            break;
          ++supergroup_iter;
        }
        // Store the new supergroup_flux
        for (unsigned int j=0; j<supergroup; ++j)
        {
          old_group_flux[i*supergroup+j] = group_flux[i*supergroup+j];
          group_flux[i*supergroup+j] = supergroup_flux[j];
        }
      }
      // Compute the convergence over all the groups
      group_conv = Compute_convergence(group_flux,old_group_flux,
          parameters.Get_n_groups());
      cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
        group_conv<<endl;
      if (group_iter>=max_group_it)
        break;
      ++group_iter;
    }

    // Apply the preconditioner to get the solution
    for (unsigned int i=0; i<parameters.Get_n_groups());
      transport_operator.Apply_preconditioner(&group_flux[i]);
    
    // Get the elapsed times
    MIP* mip(transport_operator.Get_mip());
    building_mip_time = mip->Get_building_mip_time();
    solve_mip_time = mip->Get_solve_mip_time();
    if (parameters.Get_mip_solver_type()==cg_sgs ||
        parameters.Get_mip_solver_type()==cg_ml)
      init_prec_mip_time = mip->Get_init_prec_mip_time();
  }
  else
  {
    if (parameters.Get_solver_type()!=si)
    {
      double group_conv(10.*parameters.Get_group_tol());
      Epetra_MultiVector group_flux(flux_moments_map,
          parameters.Get_n_groups());
      Epetra_MultiVector old_group_flux(group_flux);

      TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,quad[0],
          comm,flux_moments_map);
      while (group_conv>group_inner_tol)
      {
        // Loop over the supergroups.
        for (unsigned int i=0; i<parameters.Get_n_supergroups(); ++i)
        {
          unsigned int supergroup_iter(0);
          double supergroup_conv(10.*group_tol);
          Epetra_MultiVector supergroup_flux(flux_moments_map,
              parameters.Get_n_groups());
          Epetra_MultiVector old_supergroup_flux(group_flux_moments);
          while (supergroup_conv>group_tol)
          {   
            // Loop over the groups in a supergroup.
            for (unsigned int j=0; j<supergroup; ++j)
            {
              unsigned int g(i*supergroup+j);
              // Set the current group
              transport_operator.Set_group(g);
              // Compute right-hand side of GMRES (uncollided flux moments 
              // (S*inv(L)*q))
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
              solver.Iterate(parameters.Get_max_it(),inner_tol);

              // Store the new flux_moments
              old_supergroup_flux[j] = supergroup_flux[j];
              supergroup_flux[j] = flux_moments;
            }
            // Compute the convergence of the supergroup
            supergoup_conv = Compute_convergence(supergroup_flux,old_supergroup_flux,
                supergroup);
            if (parameters.Get_verbose()>0)
              cout<<"Convergence over the supergroup "<<i<<" at iteration "<<
                supergroup_iter<<": "<<supergroup_conv<<endl;
            if (supergroup_iter>=max_supergroup_it)
              break;
            ++supergroup_iter;
          }
          // Store the new supergroup_flux
          for (unsigned int j=0; j<supergroup; ++j)
          {
            old_group_flux[i*supergroup+j] = group_flux[i*supergroup+j];
            group_flux[i*supergroup+j] = supergroup_flux[j];
          }
        }
        // Compute the convergence over all the groups
        group_conv = Compute_convergence(group_flux,old_group_flux,
            parameters.Get_n_groups());
        cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
          group_conv<<endl;
        if (group_iter>=max_group_it)
          break;
        ++group_iter;
      }
      // Apply the preconditioner to get the solution
      if (parameters.Get_mip()==true)
      {
        MIP* mip(transport_operator.Get_mip());
        for (unsigned int i=0; i<parameters.Get_n_groups(); ++i)
        {
          mip->Set_group(i);
          mip->Solve(&group_flux[i]);
        }
        building_mip_time = mip->Get_building_mip_time();
        solve_mip_time = mip->Get_solve_mip_time();
        if (parameters.Get_mip_solver_type()==cg_sgs ||
            parameters.Get_mip_solver_type()==cg_ml)
          init_prec_mip_time = mip->Get_init_prec_mip_time();
      }
    }
    else
    {
      const unsigned int max_it(parameters.Get_max_it());
      double group_conv(10.*parameters.Get_group_tol());
      Epetra_MultiVector group_flux(flux_moments_map,
          parameters.Get_n_groups());
      Epetra_MultiVector old_group_flux(group_flux);

      MIP* precond(NULL);
      TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,quad[0],
          comm,flux_moments_map);

      if (parameters.Get_mip()==true)
        precond = transport_operator.Get_mip();
      while (group_conv>group_tol)
      {
        // Loop over the supergroups.
        for (unsigned int i=0; i<parameters.Get_n_supergroups(); ++i)
        {
          unsigned int supergroup_iter(0);
          double supergroup_conv(10.*group_tol);
          Epetra_MultiVector supergroup_flux(flux_moments_map,
              parameters.Get_n_groups());
          Epetra_MultiVector old_supergroup_flux(group_flux_moments);
          while (supergroup_conv>group_tol)
          {   
            // Loop over the groups in a supergroup.
            for (unsigned int j=0; j<supergroup; ++j)
            {
              unsigned int g(i*supergroup+j);
              // Set the current group
              transport_operator.Set_group(g);

              Epetra_MultiVector flux_moments_old(*flux_moments);

              for (unsigned int i=0; i<max_it; ++i)
              {
                double num(1.);
                double denom(1.);
                Epetra_BLAS blas;

                transport_operator.Compute_scattering_source(*flux_moments);
                transport_operator.Sweep(*flux_moments,true);

                if (parameters.Get_mip()==true)
                {
                  Epetra_MultiVector diff_flux(*flux_moments);
                  blas.AXPY(flux_moments_size,-1.,flux_moments_old.Values(),
                      diff_flux.Values());
                  precond->Solve(diff_flux);
                  blas.AXPY(flux_moments_size,1.,diff_flux.Values(),
                      flux_moments->Values());
                }

                Epetra_MultiVector diff_flux(*flux_moments);
                blas.AXPY(flux_moments_size,-1.,flux_moments_old.Values(),
                    diff_flux.Values());

                assert(diff_flux.Norm2(&num)==0);
                assert(flux_moments->Norm2(&denom)==0);
                convergence = num/denom;
                if (parameters.Get_verbose()>1)
                  cout<<"Convergence at iteration "<<i<<": "<<convergence<<endl;
                if (convergence<inner_tol)
                  break;

                flux_moments_old = *flux_moments;
              }

              // Store the new flux_moments
              old_supergroup_flux[j] = supergroup_flux[j];
              supergroup_flux[j] = flux_moments;
            }
            // Compute the convergence of the supergroup
            supergroup_conv = Compute_convergence(supergroup_flux,
                old_supergroup_flux,supergroup);
            if (parameters.Get_verbose()>0)
              cout<<"Convergence over the supergroup "<<i<<" at iteration "<<
                supergroup_iter<<": "<<supergroup_conv<<endl;
            if (supergroup_iter>=max_supergroup_it)
              break;
            ++supergroup_iter;
          }
          // Store the new supergroup_flux
          for (unsigned int j=0; j<supergroup; ++j)
          {
            old_group_flux[i*supergroup+j] = group_flux[i*supergroup+j];
            group_flux[i*supergroup+j] = supergroup_flux[j];
          }
        }
        // Compute the convergence over all the groups
        group_conv = Compute_convergence(group_flux,old_group_flux,
            parameters.Get_n_groups());
        cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
          group_conv<<endl;
        if (group_iter>=max_group_it)
          break;
        ++group_iter;
      }
      if (parameters.Get_mip()==true)
      {
        building_mip_time = precond->Get_building_mip_time();
        solve_mip_time = precond->Get_solve_mip_time();
        if (parameters.Get_mip_solver_type()==cg_sgs ||
            parameters.Get_mip_solver_type()==cg_ml)
          init_prec_mip_time = precond->Get_init_prec_mip_time();
      }
    }
  }
  calc_timer->stop();
  cout<<"Initialization time: "<<init_timer->totalElapsedTime()<<endl;
  cout<<"Calculation time: "<<calc_timer->totalElapsedTime()<<endl;
  if (parameters.Get_mip()==true)
  {
    cout<<"Building MIP time: "<<building_mip_time<<endl;
    if (parameters.Get_mip_solver_type()==cg_sgs || 
        parameters.Get_mip_solver_type()==cg_ml)
      cout<<"Initializing CG preconditioner time: "<<init_prec_mip_time<<endl;
    cout<<"Solving MIP time: "<<solve_mip_time<<endl;
  }
  cout<<"Total elapsed time: "<<init_timer->totalElapsedTime()+
    calc_timer->totalElapsedTime()<<endl;
}

double TRANSPORT_SOLVER::Compute_convergence(Epetra_MultiVector const &flux, 
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
  file<<quad[0]->Get_n_mom()*n_dof<<"\n";

  for (unsigned int i=0; i<=n_cells; ++i)
    file<<offset[i]<<" ";
  file<<"\n";

  const unsigned int n_points(points.size());
  for (unsigned int i=0; i<n_points; ++i)
    file<<points[i][0]<<" "<<points[i][1]<<"\n";

  const unsigned int n_values(n_dof*n_mom);
  for (unsigned int i=0; i<n_values; ++i)
    file<<(*flux_moments)[0][i]<<"\n";

  // To have the scalar flux phi_00 needs to be multiply by sqrt(weight_sum)
  const double weight_sum(parameters.Get_weight_sum());
  for (unsigned int i=0; i<n_dof; ++i)
    file<<(*flux_moments)[0][i]*sqrt(weight_sum)<<"\n";

  file.close();
}
