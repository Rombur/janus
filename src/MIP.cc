#include "MIP.hh"

MIP::MIP(unsigned int level,DOF_HANDLER* dof,PARAMETERS const* param,
    QUADRATURE const* quad,Epetra_Comm const* comm) :
  lvl(level),
  ia(NULL),
  ja(NULL),
  a(NULL),
  comm(comm),
  dof_handler(dof),
  parameters(param),
  quad(quad),
  mip_map(NULL),
  A(NULL),
  building_timer(NULL),
  init_prec_timer(NULL),
  solve_timer(NULL),
  sgs_prec(NULL),
  ml_prec(NULL)
{
  building_timer = new Teuchos::Time ("building_mip_timer");
  init_prec_timer = new Teuchos::Time ("initializing_cg_prec_timer");
  solve_timer = new Teuchos::Time ("solve_mip_timer");
}

MIP::~MIP()
{
  if (ml_prec!=NULL)
  {
    delete ml_prec;
    ml_prec = NULL;
  }

  if (sgs_prec!=NULL)
  {
    delete sgs_prec;
    sgs_prec = NULL;
  }

  if (solve_timer!=NULL)
  {
    delete solve_timer;
    solve_timer = NULL;
  }

  if (init_prec_timer!=NULL)
  {
    delete init_prec_timer;
    init_prec_timer = NULL;
  }
  
  if (building_timer!=NULL)
  {
    delete building_timer;
    building_timer = NULL;
  }

  if (A!=NULL)
  {
    delete A;
    A = NULL;
  }

  if (mip_map!=NULL)
  {
    delete mip_map;
    mip_map = NULL;
  }

  if (a!=NULL)
  {
    int n(dof_handler->Get_n_dof());
    int iprint(6);
    int ijob(-1); 
    int iter(parameters->Get_max_it());
    int nrest(1);
    double tol(parameters->Get_tolerance()/100.);
    double* f(NULL);
    double* x(NULL);
#ifdef AGMG
    dagmg_(&n,a,ja,ia,f,x,&ijob,&iprint,&nrest,&iter,&tol);
#endif 
    delete a;
    a = NULL;
  }

  if (ja!=NULL)
  {
    delete ja;
    ja = NULL;
  }

  if (ia!=NULL)
  {
    delete ia;
    ia = NULL;
  }
}

void MIP::Solve(Epetra_MultiVector &flux_moments)
{
  // Restrict the Krylov vector to build the lhs
  if (mip_map==NULL)
    mip_map = new Epetra_Map(dof_handler->Get_n_dof(),0,*comm);
  Epetra_MultiVector b(*mip_map,1);
  for (unsigned int i=0; i<dof_handler->Get_n_dof(); ++i)
    b[0][i] = flux_moments[0][i];

  // Compute the right-hand side
  Compute_rhs(flux_moments,b);

  // If A was not built yet, it is built now
  if (A==NULL && a==NULL)
  {
    int* n_entries_per_row = new int[dof_handler->Get_n_dof()];
    Compute_n_entries_per_row(n_entries_per_row);
    A = new Epetra_CrsMatrix(Copy,*mip_map,n_entries_per_row);
    Build_lhs();
  }

  // Solve the system of equation
  switch (parameters->Get_mip_solver_type())
  {
    case cg_none :
      {
        Cg_solve(flux_moments,b);
        break;
      }
    case cg_sgs :
      {
        Cg_sgs_solve(flux_moments,b);
        break;
      }
    case cg_ml :
      {
        Cg_ml_solve(flux_moments,b);
        break;
      }
#ifdef AGMG
    case agmg :
      {
        Agmg_solve(flux_moments,b);
        break;
      }
#endif
    default :
      {
        Check(false,"Unknown solver for MIP.");
      }
  }
}

void MIP::Free_ml()
{
  delete ml_prec;
  ml_prec = NULL;
}

void MIP::Compute_n_entries_per_row(int* n)
{
  vector<CELL*>::iterator cell(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
  for (; cell<cell_end; ++cell)
  {
    unsigned int coupling_dof(0);
    const unsigned int i_max((*cell)->Get_last_dof());
    const unsigned int dof_per_cell(i_max-(*cell)->Get_first_dof());
    // Loop over the edges and get the number of dof of the neighborhood
    vector<EDGE*>::iterator cell_edge((*cell)->Get_cell_edges_begin());
    vector<EDGE*>::iterator cell_edge_end((*cell)->Get_cell_edges_end());
    for (; cell_edge<cell_edge_end; ++cell_edge)
    {
      if ((*cell_edge)->Is_interior())
      {
        CELL* upwind_cell(NULL);
        if ((*cell_edge)->Get_cell_index(0)==(*cell)->Get_id())
          upwind_cell = dof_handler->Get_cell((*cell_edge)->Get_cell_index(1));
        else
          upwind_cell = dof_handler->Get_cell((*cell_edge)->Get_cell_index(0));
        coupling_dof += upwind_cell->Get_last_dof()-upwind_cell->Get_first_dof();
      }
    }
    for (unsigned int i=(*cell)->Get_first_dof(); i<i_max; ++i)
      n[i] = dof_per_cell+coupling_dof;
  }
}

void MIP::Compute_rhs(Epetra_MultiVector const &x,Epetra_MultiVector &b)
{
  Teuchos::BLAS<int,double> blas;
  Teuchos::SerialDenseMatrix<int,double> const* const D2M(quad->Get_D2M());
  vector<CELL*>::iterator cell(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
  const unsigned int n_mom(quad->Get_n_mom());

  for (; cell<cell_end; ++cell)
  {
    const unsigned int i_min((*cell)->Get_first_dof());
    const unsigned int i_max((*cell)->Get_last_dof());
    FINITE_ELEMENT const* const fe((*cell)->Get_fe());
    unsigned int dof_per_cell(fe->Get_dof_per_cell());
    Teuchos::SerialDenseVector<int,double> x_cell(dof_per_cell);
    Teuchos::SerialDenseVector<int,double> b_cell(dof_per_cell);
    Teuchos::SerialDenseMatrix<int,double> const* const mass_matrix(
        fe->Get_mass_matrix());
    vector<EDGE*>::iterator cell_edge((*cell)->Get_cell_edges_begin());
    vector<EDGE*>::iterator cell_edge_end((*cell)->Get_cell_edges_end());
    
    // Compute the volumetric terms
    for (unsigned int i=i_min; i<i_max; ++i)
      x_cell(i-i_min) = x[0][i];
    blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,
        (*cell)->Get_sigma_s(lvl,0),mass_matrix->values(),
        mass_matrix->stride(),x_cell.values(),1,0.,b_cell.values(),1);

    // Compute the reflective boundary term.
    for (; cell_edge<cell_edge_end; ++cell_edge)
    {                
      if ((*cell_edge)->Is_interior()==false)
      {
        if ((*cell_edge)->Get_bc_type()==reflective)
        {
          unsigned int edge_lid((*cell_edge)->Get_lid(0));
          const unsigned int n_dir(quad->Get_n_dir());
          Teuchos::SerialDenseVector<int,double> const* const external_normal(
              (*cell_edge)->Get_external_normal(0));
          Teuchos::SerialDenseMatrix<int,double> const* const downwind_matrix(
              fe->Get_downwind_matrix(edge_lid));
          for (unsigned int idir=0; idir<n_dir; ++idir)
          {
            Teuchos::SerialDenseVector<int,double> omega(quad->Get_omega_2d(idir));
            const double n_dot_omega(omega.dot(*external_normal));
            if (n_dot_omega>0.)
            {                                      
              const unsigned int offset(dof_handler->Get_n_dof()*n_mom+
                  dof_handler->Get_saf_map_reflective_dof((*cell)->Get_id())+
                  idir*dof_handler->Get_n_sf_per_dir());
              for (unsigned int i=0; i<dof_per_cell; ++i)
                x_cell(i) = x[0][offset+i];
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,
                  n_dot_omega*(*D2M)(0,idir),downwind_matrix->values(),
                  downwind_matrix->stride(),x_cell.values(),1,1.,b_cell.values(),1);
            }
          }
        }
      }
    }

    // Convert the result back to the Epetra_MultiVector.
    for (unsigned int i=i_min; i<i_max; ++i)
      b[0][i] = b_cell(i-i_min);
  }
}

void MIP::Build_lhs()
{
  // Start the building_timer
  building_timer->start();

  // Volumetric terms: loop over the cells
  vector<CELL*>::iterator cell(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
  for (; cell<cell_end; ++cell)
  {
    const unsigned int i_min((*cell)->Get_first_dof());
    const unsigned int i_max((*cell)->Get_last_dof());
    int* indices = new int[i_max-i_min];
    double* values = new double[i_max-i_min];
    FINITE_ELEMENT const* const fe((*cell)->Get_fe());
    Teuchos::SerialDenseMatrix<int,double> const* const mass_matrix(
        fe->Get_mass_matrix());
    Teuchos::SerialDenseMatrix<int,double> const* const stiffness_matrix(
        fe->Get_stiffness_matrix());
    for (unsigned int i=i_min; i<i_max; ++i)
      indices[i-i_min] = i;

    for (unsigned int i=i_min; i<i_max; ++i)
    {
      for (unsigned int j=i_min; j<i_max; ++j)
      {                                    
        values[j-i_min] = ((*cell)->Get_sigma_t(lvl)-(*cell)->Get_sigma_s(lvl,0))*
          (*mass_matrix)(i-i_min,j-i_min)+(*cell)->Get_diffusion_coefficient()*
          (*stiffness_matrix)(i-i_min,j-i_min);
      }

      A->InsertGlobalValues(i,i_max-i_min,&values[0],&indices[0]);
    }
    delete [] values;
    delete [] indices;
  }

  // Surfacic terms: loop over the edges
  vector<EDGE>::iterator edge(dof_handler->Get_edges_begin());
  vector<EDGE>::iterator edge_end(dof_handler->Get_edges_end());
  for (; edge<edge_end; ++edge)
  {  
    if (edge->Is_interior()==true)
    {
      const unsigned int edge_lid_0(edge->Get_lid(0));
      const unsigned int edge_lid_1(edge->Get_lid(1));
      const double normal_x(edge->Get_external_normal_component(0,0));
      const double normal_y(edge->Get_external_normal_component(0,1));
      CELL* cell(dof_handler->Get_cell(edge->Get_cell_index(0)));
      CELL* next_cell(dof_handler->Get_cell(edge->Get_cell_index(1)));
      const unsigned int cell_first_dof(cell->Get_first_dof());
      const unsigned int cell_last_dof(cell->Get_last_dof());
      const unsigned int cell_dof(cell_last_dof-cell_first_dof);
      const unsigned int next_cell_first_dof(next_cell->Get_first_dof());
      const unsigned int next_cell_last_dof(next_cell->Get_last_dof());
      const unsigned int next_cell_dof(next_cell_last_dof-next_cell_first_dof);
      FINITE_ELEMENT const* const fe(cell->Get_fe());
      FINITE_ELEMENT const* const next_fe(next_cell->Get_fe());
      // n_x dot edge_deln_matrix_m_x
      Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_m_x(
          *(fe->Get_edge_deln_matrix(edge_lid_0,0)));
      edge_deln_matrix_m_x *= normal_x;
      // n_y dot edge_deln_matrix_m_y
      Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_m_y(
          *(fe->Get_edge_deln_matrix(edge_lid_0,1)));
      edge_deln_matrix_m_y *= normal_y;
      // n dot edge_deln_matrix_m
      Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_m(
          edge_deln_matrix_m_x);
      edge_deln_matrix_m += edge_deln_matrix_m_y;
      // n_x dot edge_deln_matrix_p_x
      Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_p_x(
          *(next_fe->Get_edge_deln_matrix(edge_lid_1,0)));
      edge_deln_matrix_p_x *= normal_x;
      // n_y dot edge_deln_matrix_p_y
      Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_p_y(
          *(next_fe->Get_edge_deln_matrix(edge_lid_1,1)));
      edge_deln_matrix_p_y *= normal_y;
      // n dot edge_deln_matrix_p
      Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_p(
          edge_deln_matrix_p_x);
      edge_deln_matrix_p += edge_deln_matrix_p_y;
      // n_x dot coupling_edge_deln_matrix_m_x
      Teuchos::SerialDenseMatrix<int,double> coupling_edge_deln_matrix_m_x(
          *(fe->Get_coupling_edge_deln_matrix(edge_lid_0,0)));
      coupling_edge_deln_matrix_m_x *= normal_x;
      // n_y dot coupling_edge_deln_matrix_m_y
      Teuchos::SerialDenseMatrix<int,double> coupling_edge_deln_matrix_m_y(
          *(fe->Get_coupling_edge_deln_matrix(edge_lid_0,1)));
      coupling_edge_deln_matrix_m_y *= normal_y;
      // n dot coupling_edge_deln_matrix_m
      Teuchos::SerialDenseMatrix<int,double> coupling_edge_deln_matrix_m(
          coupling_edge_deln_matrix_m_x);
      coupling_edge_deln_matrix_m += coupling_edge_deln_matrix_m_y;
      // n_x dot coupling_edge_deln_matrix_p_x
      Teuchos::SerialDenseMatrix<int,double> coupling_edge_deln_matrix_p_x(
          *(next_fe->Get_coupling_edge_deln_matrix(edge_lid_1,0)));
      coupling_edge_deln_matrix_p_x *= normal_x;
      // n_y dot coupling_edge_deln_matrix_p_y
      Teuchos::SerialDenseMatrix<int,double> coupling_edge_deln_matrix_p_y(
          *(next_fe->Get_coupling_edge_deln_matrix(edge_lid_1,1)));
      coupling_edge_deln_matrix_p_y *= normal_y;
      // n dot coupling_edge_deln_matrix_p
      Teuchos::SerialDenseMatrix<int,double> coupling_edge_deln_matrix_p(
          coupling_edge_deln_matrix_p_x);
      coupling_edge_deln_matrix_p += coupling_edge_deln_matrix_p_y;
      const double h_m(cell->Get_orthogonal_length(edge->Get_lid(0)));
      const double h_p(next_cell->Get_orthogonal_length(edge->Get_lid(1)));
      const double D_m(cell->Get_diffusion_coefficient());
      const double D_p(next_cell->Get_diffusion_coefficient());
      const double K(Compute_penalty_coefficient(D_m,D_p,h_m,h_p,true));
      int* indices = new int[cell_last_dof-cell_first_dof];
      int* next_indices = new int[next_cell_last_dof-next_cell_first_dof];
      double* values = new double[cell_last_dof-cell_first_dof]; 
      double* next_values = new double[next_cell_last_dof-next_cell_first_dof]; 
      for (unsigned int i=cell_first_dof; i<cell_last_dof; ++i)
        indices[i-cell_first_dof] = i;
      for (unsigned int i=next_cell_first_dof; i<next_cell_last_dof; ++i)
        next_indices[i-next_cell_first_dof] = i;
      
      // First edge term ([[phi]],[[phi^*]])
      for (unsigned int i=0; i<cell_dof; ++i)
      {
        for (unsigned int j=0; j<cell_dof; ++j)
          values[j] = K*(*fe->Get_downwind_matrix(edge_lid_0))(i,j);
        A->InsertGlobalValues(i+cell_first_dof,cell_dof,&values[0],&indices[0]);
        for (unsigned int j=0; j<next_cell_dof; ++j)
          next_values[j] = -K*(*fe->Get_upwind_matrix(edge_lid_0))(i,j);
        A->InsertGlobalValues(i+cell_first_dof,next_cell_dof,&next_values[0],
            &next_indices[0]);
      }
      for (unsigned int i=0; i<next_cell_dof; ++i)
      {
        for (unsigned int j=0; j<cell_dof; ++j)
          values[j] = -K*(*next_fe->Get_upwind_matrix(edge_lid_1))(i,j);
        A->InsertGlobalValues(i+next_cell_first_dof,cell_dof,&values[0],&indices[0]);
        for (unsigned int j=0; j<next_cell_dof; ++j)
          next_values[j] = K*(*next_fe->Get_downwind_matrix(edge_lid_1))(i,j);
        A->InsertGlobalValues(i+next_cell_first_dof,next_cell_dof,&next_values[0],
            &next_indices[0]);
      }

      // Internal terms (-,-)
      for (unsigned int i=0; i<cell_dof; ++i)
      {
        for (unsigned int j=0; j<cell_dof; ++j)
          values[j] = -0.5*D_m*(edge_deln_matrix_m(i,j)+edge_deln_matrix_m(j,i));

        A->InsertGlobalValues(i+cell_first_dof,cell_dof,&values[0],&indices[0]);
      }

      // External terms (+,+)
      for (unsigned int i=0; i<next_cell_dof; ++i)
      {
        for (unsigned int j=0; j<next_cell_dof; ++j)
          next_values[j] = 0.5*D_p*(edge_deln_matrix_p(i,j)+edge_deln_matrix_p(j,i));

        A->InsertGlobalValues(i+next_cell_first_dof,next_cell_dof,&next_values[0],
            &next_indices[0]);
      }

      // Mixte terms (+,-)
      for (unsigned int i=0; i<cell_dof; ++i)
      {
        for (unsigned int j=0; j<next_cell_dof; ++j)
          next_values[j] = 0.5*(D_m*coupling_edge_deln_matrix_m(i,j)-              
              D_p*coupling_edge_deln_matrix_p(j,i));
        A->InsertGlobalValues(i+cell_first_dof,next_cell_dof,&next_values[0],
            &next_indices[0]);
      }

      // Mixte terms (-,+)
      for (unsigned int i=0; i<next_cell_dof; ++i)
      {
        for (unsigned int j=0; j<cell_dof; ++j)
          values[j] = 0.5*(-D_p*coupling_edge_deln_matrix_p(i,j)+
              D_m*coupling_edge_deln_matrix_m(j,i));

        A->InsertGlobalValues(i+next_cell_first_dof,cell_dof,&values[0],&indices[0]);
      }

      delete [] next_values;
      next_values = NULL;
      delete [] values;
      values = NULL;
      delete [] next_indices;
      next_indices = NULL;
      delete [] indices;
      indices = NULL;
    }
    else
    {
      // The edge is on the boundary. If the edge is on a reflective boundary,
      // there is nothing to do.
      if (edge->Get_bc_type()!=reflective)
      {
        const unsigned int edge_lid_0(edge->Get_lid(0));
        CELL* cell(dof_handler->Get_cell(edge->Get_cell_index(0)));
        const unsigned int cell_first_dof(cell->Get_first_dof());
        const unsigned int cell_last_dof(cell->Get_last_dof());
        const unsigned int cell_dof(cell_last_dof-cell_first_dof);
        FINITE_ELEMENT const* const fe(cell->Get_fe());
        // n_x dot edge_deln_matrix_m_x
        Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_m_x(
            *(fe->Get_edge_deln_matrix(edge_lid_0,0)));
        edge_deln_matrix_m_x *= edge->Get_external_normal_component(0,0);
        // n_y dot edge_deln_matrix_m_y
        Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_m_y(
            *(fe->Get_edge_deln_matrix(edge_lid_0,1)));
        edge_deln_matrix_m_y *= edge->Get_external_normal_component(0,1);
        // n dot edge_deln_matrix_m
        Teuchos::SerialDenseMatrix<int,double> edge_deln_matrix_m(
            edge_deln_matrix_m_x);
        edge_deln_matrix_m += edge_deln_matrix_m_y;
        const double h_m(cell->Get_orthogonal_length(edge->Get_lid(0)));
        const double D_m(cell->Get_diffusion_coefficient());
        const double K(Compute_penalty_coefficient(D_m,0.,h_m,0.,false));
        int* indices = new int[cell_last_dof-cell_first_dof];
        double* values = new double[cell_last_dof-cell_first_dof]; 
        for (unsigned int i=cell_first_dof; i<cell_last_dof; ++i)
          indices[i-cell_first_dof] = i;

        for (unsigned int i=0; i<cell_dof; ++i)
        {
          for (unsigned int j=0; j<cell_dof; ++j)
            values[j] = K*(*fe->Get_downwind_matrix(edge_lid_0))(i,j)-0.5*D_m*
              (edge_deln_matrix_m(i,j)+edge_deln_matrix_m(j,i));
          A->InsertGlobalValues(i+cell_first_dof,cell_dof,&values[0],&indices[0]);
        }
      }
    }
  }
  // Finish up, transforming the matrix entries into local numbering, to
  // optimize data transert during matrix-vector products
  A->FillComplete();
  
  // Stop the building_timer
  building_timer->stop();
}

double MIP::Compute_penalty_coefficient(double D_m,double D_p,double h_m,
    double h_p,bool interior)
{
  double k(0.);
  double k_mip(0.25);
  const double C(2.);
  // Only use first order polynomial -> c(p^+)=c(p^-)
  const double p(1.);
  const double c_p(C*p*(p+1));

  if (interior==true)
    k = c_p/2.*(D_p/h_p+D_m/h_m);
  else
    k = c_p*D_m/h_m;
  
  if (k>k_mip)
    k_mip = k;

  return k_mip;
}

void MIP::Cg_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b)
{
  Epetra_MultiVector x(*mip_map,1);
  Epetra_LinearProblem problem(A,&x,&b);
  AztecOO solver(problem);

  // No preconditioner is used
  solver.SetAztecOption(AZ_precond,AZ_none);
  // Use CG
  solver.SetAztecOption(AZ_solver,AZ_cg); 
  // Convergence criterion ||r||_2/||b||_2
  solver.SetAztecOption(AZ_conv,AZ_rhs);
  // Get the verbosity of the output
  switch (parameters->Get_verbose())
  {
    case 0 :
      {
        solver.SetAztecOption(AZ_output,AZ_none);
        break;
      }
    case 1 :
      {
        solver.SetAztecOption(AZ_output,AZ_summary);
        break;
      }
    case 2 :
      {
        solver.SetAztecOption(AZ_output,AZ_all);
        break;
      }
    default :
      {
        solver.SetAztecOption(AZ_output,AZ_none);
      }
  }

  // Start solve_timer
  solve_timer->start();
  // Solve the MIP equation
  solver.Iterate(parameters->Get_max_it(),parameters->Get_tolerance()/100.);
  // Stop solve_timer
  solve_timer->stop();

  // Project the solution
  Project_solution(flux_moments,x);
}

void MIP::Cg_sgs_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b)
{
  Epetra_MultiVector x(*mip_map,1);
  Epetra_LinearProblem problem(A,&x,&b);
  AztecOO solver(problem);
  
  if (sgs_prec==NULL)
  {
    // Start init_prec_timer
    init_prec_timer->start();
    Teuchos::ParameterList sgs_list;
    sgs_list.set("relaxation: type", "symmetric Gauss-Seidel");
    sgs_prec = new Ifpack_PointRelaxation(A);
    sgs_prec->SetParameters(sgs_list);
    sgs_prec->Initialize();
    sgs_prec->Compute();
    // Stop init_prec_timer
    init_prec_timer->stop();
  }
  // Set the preconditioner to symmetric Gauss-Seidel
  solver.SetPrecOperator(sgs_prec);
  // Use CG
  solver.SetAztecOption(AZ_solver,AZ_cg); 
  // Convergence criterion ||r||_2/||b||_2
  solver.SetAztecOption(AZ_conv,AZ_rhs);
  // Get the verbosity of the output
  switch (parameters->Get_verbose())
  {
    case 0 :
      {
        solver.SetAztecOption(AZ_output,AZ_none);
        break;
      }
    case 1 :
      {
        solver.SetAztecOption(AZ_output,AZ_summary);
        break;
      }
    case 2 :
      {
        solver.SetAztecOption(AZ_output,AZ_all);
        break;
      }
    default :
      {
        solver.SetAztecOption(AZ_output,AZ_none);
      }
  }

  // Start solve_timer
  solve_timer->start();
  // Solve the MIP equation
  solver.Iterate(parameters->Get_max_it(),parameters->Get_tolerance()/100.);
  // Stop solve_timer
  solve_timer->stop();

  // Project the solution
  Project_solution(flux_moments,x);
}

void MIP::Cg_ml_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b)
{
  Epetra_MultiVector x(*mip_map,1);
  Epetra_LinearProblem problem(A,&x,&b);
  AztecOO solver(problem);

  if (ml_prec==NULL)
  {
    // Start init_prec_timer
    init_prec_timer->start();
    Teuchos::ParameterList ml_list;
    // Set default values for classical smoothed aggregation in ml_list
    ML_Epetra::SetDefaults("SA",ml_list);

    // Overwrite some of the parameters
    if (parameters->Get_aggregation_type()==uncoupled)
      ml_list.set("aggregation: type","Uncoupled");
    else 
      if (parameters->Get_aggregation_type()==mis) 
        ml_list.set("aggregation: type","MIS");
    if (parameters->Get_verbose()==2)
      ml_list.set("ML output",10);

    // Construct the ML object and compute the hierarchy
    ml_prec = new ML_Epetra::MultiLevelPreconditioner(*A,ml_list,true);
    // Stop init_prec_timer
    init_prec_timer->stop();
  }
  // Set the AMG as preconditioner
  solver.SetPrecOperator(ml_prec);
  // Use CG
  solver.SetAztecOption(AZ_solver,AZ_cg); 
  // Convergence criterion ||r||_2/||b||_2
  solver.SetAztecOption(AZ_conv,AZ_rhs);
  // Get the verbosity of the output
  switch (parameters->Get_verbose())
  {
    case 0 :
      {
        solver.SetAztecOption(AZ_output,AZ_none);
        break;
      }
    case 1 :
      {
        solver.SetAztecOption(AZ_output,AZ_summary);
        break;
      }
    case 2 :
      {
        solver.SetAztecOption(AZ_output,AZ_all);
        break;
      }
    default :
      {
        solver.SetAztecOption(AZ_output,AZ_none);
      }
  }

  // Start solve_timer
  solve_timer->start();
  // Solve the MIP equation
  solver.Iterate(parameters->Get_max_it(),parameters->Get_tolerance()/100.);
  // Stop solve_timer
  solve_timer->stop();

  // Project the solution
  Project_solution(flux_moments,x);
}

#ifdef AGMG
void MIP::Agmg_solve(Epetra_MultiVector &flux_moments,Epetra_MultiVector &b)
{
  int n(dof_handler->Get_n_dof());
  int iprint(6);
  int ijob(2); 
  int iter(parameters->Get_max_it());
  int nrest(1);
  double tol(parameters->Get_tolerance()/100.);
  double* f;
  double* x;
  Epetra_MultiVector epetra_x(*mip_map,1);
  x = new double [n];
  f = new double [n];
  
  // Convert the right-hand side to fortran
  Convert_rhs_to_fortran(b,f,n);

  // If A has not been converted to fortran yet, we do it now
  if (a==NULL)
  {
    int ijob_1(1);
    Convert_lhs_to_fortran(ia,ja,a,n,A->NumGlobalNonzeros());
    dagmg_(&n,a,ja,ia,f,x,&ijob_1,&iprint,&nrest,&iter,&tol);
    // Stop init_prec_timer
    init_prec_timer->stop();
  }

  // Start solve_timer
  solve_timer->start();
  // Solve the MIP equation
  dagmg_(&n,a,ja,ia,f,x,&ijob,&iprint,&nrest,&iter,&tol);
  // Stop solve_timer
  solve_timer->stop();

  // Convert the righ-hand side to Epetra
  Convert_rhs_to_epetra(epetra_x,x,n);
  
  // Project the solution
  Project_solution(flux_moments,epetra_x);

  delete [] f;
  delete [] x;
}
#endif

void MIP::Project_solution(Epetra_MultiVector &flux_moments,
    Epetra_MultiVector const &x)
{                    
  const unsigned int i_max(dof_handler->Get_n_dof());
  const unsigned int n_dir(quad->Get_n_dir());
  const unsigned int n_mom(quad->Get_n_mom());
  Teuchos::SerialDenseMatrix<int,double> const* const M2D(quad->Get_M2D());

  if (parameters->Get_solver_type()==si)
  {
    for (unsigned int i=0; i<i_max; ++i)
      flux_moments[0][i] = x[0][i];
    // Add the correction on the angular flux when there is reflective
    // boundary
    vector<EDGE>::iterator edge(dof_handler->Get_edges_begin());
    vector<EDGE>::iterator edge_end(dof_handler->Get_edges_end());
    for (; edge<edge_end; ++edge)
    {  
      if (edge->Is_interior()==false)
      {
        if (edge->Get_bc_type()==reflective)
        {
          CELL* cell(dof_handler->Get_cell(edge->Get_cell_index(0)));
          const unsigned int dof_per_cell(cell->Get_last_dof()-
              cell->Get_first_dof());
          const double D(cell->Get_diffusion_coefficient());
          for (unsigned int idir=0; idir<n_dir; ++idir)
          {
            const unsigned int offset_1(dof_handler->Get_n_dof()*n_mom+
                dof_handler->Get_saf_map_reflective_dof(cell->Get_id())+
                idir*dof_handler->Get_n_sf_per_dir());
            const unsigned int offset_2(cell->Get_first_dof());
            const double M2D0((*M2D)(idir,0));
            if (n_mom==1)
            {     
              for (unsigned int j=0; j<dof_per_cell; ++j)
                flux_moments[0][offset_1+j] = M2D0*x[0][offset_2+j];
            }
            else
            {
              const unsigned int edge_lid(edge->Get_lid(0));
              Teuchos::SerialDenseVector<int,double> omega(quad->Get_omega_2d(idir));
              const double omega_0(omega(0));
              const double omega_1(omega(1));
              const double M2D1((*M2D)(idir,1));
              const double M2D2((*M2D)(idir,2));
              FINITE_ELEMENT const* const fe(cell->Get_fe());
              Teuchos::BLAS<int,double> blas;
              Teuchos::SerialDenseMatrix<int,double> const* const x_grad(
                  fe->Get_x_grad_matrix(edge_lid));
              Teuchos::SerialDenseMatrix<int,double> const* const y_grad(
                  fe->Get_y_grad_matrix(edge_lid));
              Teuchos::SerialDenseVector<int,double> x_grad_correc(dof_per_cell);
              Teuchos::SerialDenseVector<int,double> y_grad_correc(dof_per_cell);
              Teuchos::SerialDenseVector<int,double> flux_correc(dof_per_cell);
              for (unsigned int j=0; j<dof_per_cell; ++j)
                flux_correc.values()[j] = x[0][offset_2+j];
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,1.,
                  x_grad->values(),x_grad->stride(),flux_correc.values(),1,0.,
                  x_grad_correc.values(),1);
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,1.,
                  y_grad->values(),y_grad->stride(),flux_correc.values(),1,0.,
                  y_grad_correc.values(),1);
              for (unsigned int j=0; j<dof_per_cell; ++j)
                flux_moments[0][offset_1+j] = M2D0*x[0][offset_2+j]-
                  M2D1*D*omega_0*x_grad_correc(j)-M2D2*D*omega_1*y_grad_correc(j);
            }
          }
        }
      }
    }
  }           
  else
  {
    for (unsigned int i=0; i<i_max; ++i)
      flux_moments[0][i] += x[0][i];
    // Add the correction on the angular flux when there is reflective
    // boundary
    vector<EDGE>::iterator edge(dof_handler->Get_edges_begin());
    vector<EDGE>::iterator edge_end(dof_handler->Get_edges_end());
    for (; edge<edge_end; ++edge)
    {  
      if (edge->Is_interior()==false)
      {
        if (edge->Get_bc_type()==reflective)
        {
          CELL* cell(dof_handler->Get_cell(edge->Get_cell_index(0)));
          const unsigned int dof_per_cell(cell->Get_last_dof()-
              cell->Get_first_dof());
          const double D(cell->Get_diffusion_coefficient());
          for (unsigned int idir=0; idir<n_dir; ++idir)
          {
            const unsigned int offset_1(dof_handler->Get_n_dof()*n_mom+
                dof_handler->Get_saf_map_reflective_dof(cell->Get_id())+
                idir*dof_handler->Get_n_sf_per_dir());
            const unsigned int offset_2(cell->Get_first_dof());
            const double M2D0((*M2D)(idir,0));
            if (n_mom==1)
            {     
              for (unsigned int j=0; j<dof_per_cell; ++j)
                flux_moments[0][offset_1+j] += M2D0*x[0][offset_2+j];
            }
            else
            {
              const unsigned int edge_lid(edge->Get_lid(0));
              Teuchos::SerialDenseVector<int,double> omega(quad->Get_omega_2d(idir));
              const double omega_0(omega(0));
              const double omega_1(omega(1));
              const double M2D1((*M2D)(idir,1));
              const double M2D2((*M2D)(idir,2));
              FINITE_ELEMENT const* const fe(cell->Get_fe());
              Teuchos::BLAS<int,double> blas;
              Teuchos::SerialDenseMatrix<int,double> const* const x_grad(
                  fe->Get_x_grad_matrix(edge_lid));
              Teuchos::SerialDenseMatrix<int,double> const* const y_grad(
                  fe->Get_y_grad_matrix(edge_lid));
              Teuchos::SerialDenseVector<int,double> x_grad_correc(dof_per_cell);
              Teuchos::SerialDenseVector<int,double> y_grad_correc(dof_per_cell);
              Teuchos::SerialDenseVector<int,double> flux_correc(dof_per_cell);
              for (unsigned int j=0; j<dof_per_cell; ++j)
                flux_correc.values()[j] = x[0][offset_2+j];
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,1.,
                  x_grad->values(),x_grad->stride(),flux_correc.values(),1,0.,
                  x_grad_correc.values(),1);
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,1.,
                  y_grad->values(),y_grad->stride(),flux_correc.values(),1,0.,
                  y_grad_correc.values(),1);
              for (unsigned int j=0; j<dof_per_cell; ++j)
                flux_moments[0][offset_1+j] += M2D0*x[0][offset_2+j]-
                  M2D1*D*omega_0*x_grad_correc(j)-M2D2*D*omega_1*y_grad_correc(j);
            }
          }
        }
      }
    }
  }
}
  
void MIP::Convert_lhs_to_fortran(int* &ia,int* &ja,double* &a,unsigned int n_dof,
    unsigned int nnz)
{
  // Fortran array starts at 1 not 0
  int offset(1);
  a = new double[nnz];
  ja = new int[nnz];
  ia = new int[n_dof+1];
  ia[0] = offset;

  for (unsigned int i=0; i<n_dof; ++i)
  {
    // Copy A to a, ia and ja
    int nnz_row(A->NumGlobalEntries(i));
    double* val = new double [nnz_row];
    int* ind = new int [nnz_row];
    int n_entries;
    
    A->ExtractGlobalRowCopy(i,nnz_row,n_entries,val,ind);

    for (int j=0; j<n_entries; ++j)
    {
      a[offset-1+j] = val[j];
      // Fortran array starts at 1 not 0
      ja[offset-1+j] = ind[j]+1;
    }

    offset += n_entries;
    ia[i+1] = offset;

    delete [] ind;
    delete [] val;
  }

  delete A;
  A = NULL;
}

void MIP::Convert_rhs_to_fortran(Epetra_MultiVector const &b,double* f,
    const unsigned int n_dof) const
{
  // Copy b to f
  for (unsigned int i=0; i<n_dof; ++i)
    f[i] = b[0][i];
}

void MIP::Convert_rhs_to_epetra(Epetra_MultiVector &b,double* f,unsigned int n_dof) 
  const
{
  // Copy f to b
  for (unsigned int i=0; i<n_dof; ++i)
    b[0][i] = f[i];
}
