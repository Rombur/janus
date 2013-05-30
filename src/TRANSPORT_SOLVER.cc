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

#include "TRANSPORT_SOLVER.hh"

TRANSPORT_SOLVER::TRANSPORT_SOLVER(string* g_inputfile,string* p_inputfile,
    string* xs_inputfile,string* outputfile,Epetra_MpiComm* mpi_comm) :
  init_timer(NULL),
  calc_timer(NULL),
  outputfile(outputfile),
  comm(mpi_comm),
  group_flux(NULL),
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
  parameters.Read_parameters(triangulation.Get_n_sources());
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
            parameters.Get_permutation_type(),false);
        break;
      }
    case regular_exs :
      {
        cross_sections.Read_regular_cross_sections(triangulation.Get_n_materials(),
            parameters.Get_permutation_type(),true);
        break;
      }
    case cepxs :
      {         
        cross_sections.Read_cepxs_cross_sections(triangulation.Get_n_materials(),
            parameters.Get_permutation_type());
        break;
      }
    default :
      {
        Check(false,"Unknown type of cross section file.");
      }
  }
  cross_sections.Apply_ang_lvls_and_tc(parameters.Get_multigrid(),
      parameters.Get_transport_correction(),parameters.Get_optimal_tc(),
      triangulation.Get_n_materials(),parameters.Get_sn_order());
  cout<<"Done reading the cross sections"<<endl;

  // Build the quadratures
  cout<<"Start building the quadratures"<<endl;
  unsigned int n_lvl(parameters.Get_n_levels());
  unsigned int tmp_sn(parameters.Get_sn_order());
  unsigned int tmp_L_max(cross_sections.Get_L_max());
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


  // Instantiate the flux moments map and group_flux
  flux_moments_size = dof_handler->Get_n_dof()*quad[0]->Get_n_mom()+
    dof_handler->Get_n_sf_per_dir()*quad[0]->Get_n_dir();
  flux_moments_map = new Epetra_Map(flux_moments_size,0,*comm);
  group_flux = new Epetra_MultiVector(*flux_moments_map,
      cross_sections.Get_n_groups());
  init_timer->stop();
}

TRANSPORT_SOLVER::~TRANSPORT_SOLVER()
{
  if (group_flux!=NULL)
  {
    delete group_flux;
    group_flux = NULL;
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
  const unsigned int max_group_it(parameters.Get_max_group_it());
  const unsigned int n_groups(cross_sections.Get_n_groups());
  const unsigned int max_supergroup_it(parameters.Get_max_supergroup_it());
  const unsigned int supergroup(cross_sections.Get_n_grps_in_supergrp());
  double building_mip_time(0.);
  double init_prec_mip_time(0.);
  double solve_mip_time(0.);
  const double inner_tol(parameters.Get_inner_tolerance());
  const double group_tol(parameters.Get_group_tolerance());
  if (parameters.Get_multigrid()==true)
  {
    const unsigned int lvl(0);
    const unsigned int max_lvl(parameters.Get_n_levels()-1);
    double group_conv(10.*group_tol);
    Epetra_MultiVector old_group_flux(*group_flux);
    Epetra_MultiVector supergroup_flux(*flux_moments_map,supergroup);
    Epetra_MultiVector old_supergroup_flux(supergroup_flux);
    Epetra_MultiVector flux_moments(*flux_moments_map,1);

    TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,&quad,comm,
        flux_moments_map,lvl,max_lvl,n_groups);
    while (group_conv>parameters.Get_group_tolerance())
    {
      // Loop over the supergroups.
      for (unsigned int i=0; i<cross_sections.Get_n_supergroups(); ++i)
      {
        unsigned int supergroup_iter(0);
        double supergroup_conv(10.*group_tol);
        while (supergroup_conv>group_tol)
        {   
          // Loop over the groups in a supergroup.
          for (unsigned int j=0; j<supergroup; ++j)
          {
            unsigned int g(i*supergroup+j);
            // Set the current group
            transport_operator.Set_group(g);
            // Compute right-hand side of GMRES (uncollided flux moments (S*inv(L)*q))
            Epetra_MultiVector rhs(flux_moments);
            transport_operator.Sweep(rhs,group_flux);

            Epetra_LinearProblem problem(&transport_operator,&flux_moments,&rhs);

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
            solver.Iterate(parameters.Get_max_inner_it(),inner_tol);

            // Store the new flux_moments
            copy(supergroup_flux[j],supergroup_flux[j]+flux_moments_size,
                old_supergroup_flux[j]);
            copy(flux_moments[0],flux_moments[0]+flux_moments_size,
                supergroup_flux[j]);
          }
          // Compute the convergence of the supergroup
          supergroup_conv = Compute_convergence(supergroup_flux,
              old_supergroup_flux,supergroup);
          if (parameters.Get_verbose()>0)
          cout<<"Convergence over the supergroup "<<i<<" at iteration "<<
            supergroup_iter<<": "<<supergroup_conv<<endl;
          ++supergroup_iter;
          if (supergroup_iter>=max_supergroup_it)
            break;
        }
        // Store the new supergroup_flux
        for (unsigned int j=0; j<supergroup; ++j)
        {
          copy((*group_flux)[i*supergroup+j],(*group_flux)[i*supergroup+j]+
              flux_moments_size,old_group_flux[i*supergroup+j]);
          copy(supergroup_flux[i*supergroup+j],supergroup_flux[i*supergroup+j]+
              flux_moments_size,(*group_flux)[i*supergroup+j]);
        }
      }
      // Compute the convergence over all the groups
      group_conv = Compute_convergence(*group_flux,old_group_flux,n_groups);
      cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
        group_conv<<endl;
      ++group_iter;
      if (group_iter>=max_group_it)
        break;
    }

    // Apply the preconditioner to get the solution
    for (unsigned int g=0; g<n_groups; ++g)
    {
      Epetra_MultiVector tmp(*flux_moments_map,1);
      tmp[0] = (*group_flux)[g];
      transport_operator.Set_group(g);
      transport_operator.Apply_preconditioner(tmp);
      copy(tmp[0],tmp[0]+flux_moments_size,(*group_flux)[g]);
    }
    
    // Get the elapsed times
    MIP* mip(transport_operator.Get_mip());
    building_mip_time = mip->Get_building_mip_time();
    solve_mip_time = mip->Get_solve_mip_time();
    if (parameters.Get_mip_solver_type()==cg_ssor ||
        parameters.Get_mip_solver_type()==cg_ml)
      init_prec_mip_time = mip->Get_init_prec_mip_time();
  }
  else
  {
    if (parameters.Get_solver_type()!=si)
    {
      double group_conv(10.*group_tol);
      Epetra_MultiVector old_group_flux(*group_flux);
      Epetra_MultiVector supergroup_flux(*flux_moments_map,n_groups);
      Epetra_MultiVector old_supergroup_flux(supergroup_flux);
      Epetra_MultiVector flux_moments(*flux_moments_map,1);

      TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,quad[0],
          comm,flux_moments_map,n_groups);
      while (group_conv>group_tol)
      {
        // Loop over the supergroups.
        for (unsigned int i=0; i<cross_sections.Get_n_supergroups(); ++i)
        {
          unsigned int supergroup_iter(0);
          double supergroup_conv(10.*group_tol);
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
              Epetra_MultiVector rhs(flux_moments);
              transport_operator.Sweep(rhs,group_flux);

              Epetra_LinearProblem problem(&transport_operator,&flux_moments,&rhs);

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
              solver.Iterate(parameters.Get_max_inner_it(),inner_tol);

              // Store the new flux_moments
              copy(supergroup_flux[j],supergroup_flux[j]+flux_moments_size,
                  old_supergroup_flux[j]);
              copy(flux_moments[0],flux_moments[0]+flux_moments_size,
                  supergroup_flux[j]);
            }
            // Compute the convergence of the supergroup
            supergroup_conv = Compute_convergence(supergroup_flux,old_supergroup_flux,
                supergroup);
            if (parameters.Get_verbose()>0)
              cout<<"Convergence over the supergroup "<<i<<" at iteration "<<
                supergroup_iter<<": "<<supergroup_conv<<endl;
            ++supergroup_iter;
            if (supergroup_iter>=max_supergroup_it)
              break;
          }
          // Store the new supergroup_flux
          for (unsigned int j=0; j<supergroup; ++j)
          {
            copy((*group_flux)[i*supergroup+j],(*group_flux)[i*supergroup+j]+
                flux_moments_size,old_group_flux[i*supergroup+j]);
            copy(supergroup_flux[j],supergroup_flux[j]+flux_moments_size,
                (*group_flux)[i*supergroup+j]);
          }
        }
        // Compute the convergence over all the groups
        group_conv = Compute_convergence(*group_flux,old_group_flux,n_groups);
        cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
          group_conv<<endl;
        ++group_iter;
        if (group_iter>=max_group_it)
          break;
      }
      // Apply the preconditioner to get the solution
      if (parameters.Get_mip()==true)
      {
        MIP* mip(transport_operator.Get_mip());
        for (unsigned int g=0; g<n_groups; ++g)
        {
          Epetra_MultiVector tmp(*flux_moments_map,1);
          tmp[0] = (*group_flux)[g];
          mip->Set_group(g);
          mip->Solve(tmp);
          copy(tmp[0],tmp[0]+flux_moments_size,(*group_flux)[g]);
        }
        building_mip_time = mip->Get_building_mip_time();
        solve_mip_time = mip->Get_solve_mip_time();
        if (parameters.Get_mip_solver_type()==cg_ssor ||
            parameters.Get_mip_solver_type()==cg_ml)
          init_prec_mip_time = mip->Get_init_prec_mip_time();
      }
    }
    else
    {
      const unsigned int max_it(parameters.Get_max_inner_it());
      double group_conv(10.*group_tol);
      Epetra_MultiVector old_group_flux(*group_flux);
      Epetra_MultiVector supergroup_flux(*flux_moments_map,n_groups);
      Epetra_MultiVector old_supergroup_flux(supergroup_flux);
      Epetra_MultiVector flux_moments(*flux_moments_map,1);

      MIP* precond(NULL);
      TRANSPORT_OPERATOR transport_operator(dof_handler,&parameters,quad[0],
          comm,flux_moments_map,n_groups);

      if (parameters.Get_mip()==true)
        precond = transport_operator.Get_mip();
      while (group_conv>group_tol)
      {
        // Loop over the supergroups.
        for (unsigned int i=0; i<cross_sections.Get_n_supergroups(); ++i)
        {
          unsigned int supergroup_iter(0);
          double supergroup_conv(10.*group_tol);
          while (supergroup_conv>group_tol)
          {   
            // Loop over the groups in a supergroup.
            for (unsigned int j=0; j<supergroup; ++j)
            {
              unsigned int g(i*supergroup+j);
              // Set the current group
              transport_operator.Set_group(g);

              Epetra_MultiVector flux_moments_old(flux_moments);

              for (unsigned int i=0; i<max_it; ++i)
              {
                double num(1.);
                double denom(1.);
                Epetra_BLAS blas;

                transport_operator.Compute_scattering_source(flux_moments);
                Epetra_MultiVector rhs(flux_moments);
                transport_operator.Sweep(rhs,group_flux);

                if (parameters.Get_mip()==true)
                {
                  Epetra_MultiVector diff_flux(flux_moments);
                  blas.AXPY(flux_moments_size,-1.,flux_moments_old.Values(),
                      diff_flux.Values());
                  precond->Solve(diff_flux);
                  blas.AXPY(flux_moments_size,1.,diff_flux.Values(),
                      flux_moments.Values());
                }

                Epetra_MultiVector diff_flux(flux_moments);
                blas.AXPY(flux_moments_size,-1.,flux_moments_old.Values(),
                    diff_flux.Values());

                Check(diff_flux.Norm2(&num)!=0,string ("Difference of flux between two SI iterations is 0."));
                Check(flux_moments.Norm2(&denom)!=0,string ("Flux moments after the SI iteration is 0."));
                double convergence(num/denom);
                if (parameters.Get_verbose()>1)
                  cout<<"Convergence at iteration "<<i<<": "<<convergence<<endl;
                if (convergence<inner_tol)
                  break;

                flux_moments_old = flux_moments;
              }

              // Store the new flux_moments
              copy(supergroup_flux[j],supergroup_flux[j]+flux_moments_size,
                  old_supergroup_flux[j]);
              copy(flux_moments[0],flux_moments[0]+flux_moments_size,
                  supergroup_flux[j]);
            }
            // Compute the convergence of the supergroup
            supergroup_conv = Compute_convergence(supergroup_flux,
                old_supergroup_flux,supergroup);
            if (parameters.Get_verbose()>0)
              cout<<"Convergence over the supergroup "<<i<<" at iteration "<<
                supergroup_iter<<": "<<supergroup_conv<<endl;
            ++supergroup_iter;
            if (supergroup_iter>=max_supergroup_it)
              break;
          }
          // Store the new supergroup_flux
          for (unsigned int j=0; j<supergroup; ++j)
          {
            copy((*group_flux)[i*supergroup+j],(*group_flux)[i*supergroup+j]+
                flux_moments_size,old_group_flux[i*supergroup+j]);
            copy(supergroup_flux[j],supergroup_flux[j]+flux_moments_size,
                (*group_flux)[i*supergroup+j]);
          }
        }
        // Compute the convergence over all the groups
        group_conv = Compute_convergence(*group_flux,old_group_flux,n_groups);
        cout<<"Convergence over all groups at iteration "<<group_iter<<": "<<
          group_conv<<endl;
        ++group_iter;
        if (group_iter>=max_group_it)
          break;
      }
      if (parameters.Get_mip()==true)
      {
        building_mip_time = precond->Get_building_mip_time();
        solve_mip_time = precond->Get_solve_mip_time();
        if (parameters.Get_mip_solver_type()==cg_ssor ||
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
    if (parameters.Get_mip_solver_type()==cg_ssor || 
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

void TRANSPORT_SOLVER::Write_in_file() const
{
  const unsigned int n_cells(dof_handler->Get_n_cells());
  const unsigned int n_dof(dof_handler->Get_n_dof());
  const unsigned int n_groups(cross_sections.Get_n_groups());
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
  file<<n_groups<<" ";
  file<<quad[0]->Get_n_mom()*n_dof<<" ";

  for (unsigned int i=0; i<=n_cells; ++i)
    file<<offset[i]<<" ";
  file<<"\n";

  const unsigned int n_points(points.size());
  for (unsigned int i=0; i<n_points; ++i)
    file<<points[i][0]<<" "<<points[i][1]<<"\n";

  const unsigned int n_values(n_dof*n_mom);
  const double weight_sum(parameters.Get_weight_sum());
  // Loop over the group
  for (unsigned int g=0; g<n_groups; ++g)
    for (unsigned int i=0; i<n_values; ++i)
      file<<(*group_flux)[g][i]<<"\n";

  // To have the scalar flux phi_00 needs to be multiply by sqrt(weight_sum)
  for (unsigned int g=0; g<n_groups; ++g)
    for (unsigned int i=0; i<n_dof; ++i)
      file<<(*group_flux)[g][i]*sqrt(weight_sum)<<"\n";

  // If sigma_e is given, the dose is output
  if (cross_sections.Sigma_e_exist()==false)
    file<<"false\n";
  else
  {
    file<<"true\n";
    for (cell_it = dof_handler->Get_mesh_begin();cell_it<cell_end; ++cell_it)
    {
      for (unsigned int i=0; i<n_values; ++i)
      {
        double dose(0.);
        for (unsigned int g=0; g<n_groups; ++g)
          dose += (*group_flux)[g][i]*(*cell_it)->Get_sigma_e(g);
        file<<dose<<"\n";
      }
    }
  }


  file.close();
}
