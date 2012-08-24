#include "TRANSPORT_OPERATOR.hh"

TRANSPORT_OPERATOR::TRANSPORT_OPERATOR(DOF_HANDLER* dof,
    PARAMETERS const* param, QUADRATURE* quad,Epetra_Comm const* comm,
    Epetra_Map const* flux_moments_map) :
  Epetra_Operator(),
  lvl(0),
  max_lvl(0),
  n_dof(dof->Get_n_dof()),
  comm(comm),
  flux_moments_map(flux_moments_map),
  saf_map(NULL),
  saf(NULL),
  scattering_src(NULL),
  dof_handler(dof),
  precond(NULL),
  param(param),
  quad(quad)
{
  if (dof_handler->Get_n_sf_per_dir()>0)
  {
    saf_map = new Epetra_Map(dof_handler->Get_n_sf_per_dir()*quad->Get_n_dir(),
        0,*comm);
    saf = new Epetra_MultiVector(*saf_map,1);
  }

  // Adapt the size of teuchos_vector
  teuchos_vector.resize(dof_handler->Get_max_dof_per_cell());
  for (unsigned int i=0; i<=dof_handler->Get_max_dof_per_cell(); ++i)
    teuchos_vector[i] = new Teuchos::SerialDenseVector<int,double> (i);

  scattering_src = new vector<Teuchos::SerialDenseVector<int,double> >
    (quad->Get_n_mom(),Teuchos::SerialDenseVector<int,double> 
     (dof_handler->Get_n_dof()));
  if (param->Get_mip()==true)
  {
    if (param->Get_transport_correction()==true)
      precond = new MIP (1,dof,param,quad,comm);
    else
      precond = new MIP (0,dof,param,quad,comm);
  }
}

TRANSPORT_OPERATOR::TRANSPORT_OPERATOR(DOF_HANDLER* dof,
    PARAMETERS const* param,vector<QUADRATURE*> const* quad_vector,
    Epetra_Comm const* comm,Epetra_Map const* flux_moments_map,unsigned int level,
    unsigned int max_level,MIP* preconditioner) :
  Epetra_Operator(),
  lvl(level),
  max_lvl(max_level),
  n_dof(dof->Get_n_dof()),
  comm(comm),
  flux_moments_map(flux_moments_map),
  saf_map(NULL),
  saf(NULL),
  scattering_src(NULL),
  dof_handler(dof),
  precond(preconditioner),
  param(param),
  quad((*quad_vector)[lvl]),
  quad_vector(quad_vector)
{
  // Adapt the size of teuchos_vector
  teuchos_vector.resize(dof_handler->Get_max_dof_per_cell());
  for (unsigned int i=0; i<=dof_handler->Get_max_dof_per_cell(); ++i)
    teuchos_vector[i] = new Teuchos::SerialDenseVector<int,double> (i);

  scattering_src = new vector<Teuchos::SerialDenseVector<int,double> >
    (quad->Get_n_mom(),Teuchos::SerialDenseVector<int,double> 
     (dof_handler->Get_n_dof()));
  if (lvl==0)
    precond = new MIP (max_lvl,dof,param,quad,comm);
}

TRANSPORT_OPERATOR::~TRANSPORT_OPERATOR()
{
  if (precond!=NULL && lvl==0)
  {
    delete precond;
    precond = NULL;
  }

  if (scattering_src!=NULL)
  {
    delete scattering_src;
    scattering_src = NULL;
  }

  for (unsigned int i=0; i<=dof_handler->Get_max_dof_per_cell(); ++i)
  {
    delete teuchos_vector[i];
    teuchos_vector[i] = NULL;
  }

  if (saf!=NULL)
  {
    delete saf;
    saf = NULL;
  }

  if (saf_map!=NULL)
  {
    delete saf_map;
    saf_map = NULL;
  }
}

int TRANSPORT_OPERATOR::Apply(Epetra_MultiVector const &x,Epetra_MultiVector &y) const
{
  y = x;

  if (param->Get_multigrid()==true)
  {
    if (lvl!=max_lvl-1)
    {
      Epetra_MultiVector z(y);
      Epetra_Map coarse_map(n_dof*(*quad_vector)[lvl+1]->Get_n_mom(),0,*comm);
      TRANSPORT_OPERATOR coarse_transport(dof_handler,param,quad_vector,comm,
          &coarse_map,lvl+1,max_lvl,precond);

      Epetra_MultiVector coarse_y(coarse_transport.Restrict_vector(y));
      Epetra_MultiVector coarse_x(coarse_y);
      coarse_transport.Apply(coarse_x,coarse_y);

      // Project y on z. 
      Project_vector(z,coarse_y);
      y = z;

      // Compute the scattering source
      Compute_scattering_source(y);
      Sweep(y);      
      if (lvl==0)
      {
        for (unsigned int i=0; i<n_dof*quad->Get_n_mom(); ++i)
          y[0][i] = z[0][i]-y[0][i];
      }
    }
    else
    {
      // Apply MIP
      precond->Solve(y);

      // Compute the scattering source
      Compute_scattering_source(y);
      Sweep(y);
    }
  }
  else
  {
    if (param->Get_mip()==true)
      precond->Solve(y);
    Epetra_MultiVector z(y);
    
    // Compute the scattering source
    Compute_scattering_source(y);
    Sweep(y);

    for (int i=0; i<y.MyLength(); ++i)
      y[0][i] = z[0][i]-y[0][i];
  }  

  return 0;
}

void TRANSPORT_OPERATOR::Apply_preconditioner(Epetra_MultiVector &x) 
{
  if (lvl!=max_lvl-1)
  {
    Epetra_MultiVector y(x);
    Epetra_Map coarse_map(n_dof*(*quad_vector)[lvl+1]->Get_n_mom(),0,*comm);
    TRANSPORT_OPERATOR coarse_transport(dof_handler,param,quad_vector,comm,
        &coarse_map,lvl+1,max_lvl,precond);

    Epetra_MultiVector coarse_y(coarse_transport.Restrict_vector(y));
    coarse_transport.Apply_preconditioner(coarse_y);

    // Project y on x. 
    Project_vector(x,coarse_y);

    if (lvl!=0)
    {
      // Compute the scattering source
      Compute_scattering_source(x);
      Sweep(x);      
    }
  }
  else
  {
    // Apply MIP
    precond->Solve(x);

    // Compute the scattering source
    Compute_scattering_source(x);
    Sweep(x);
  }
}

void TRANSPORT_OPERATOR::Compute_scattering_source(Epetra_MultiVector const &x) const
{
  const unsigned int n_mom(quad->Get_n_mom());
  // Reinitialize the scattering source
  for (unsigned int i=0; i<n_mom; ++i)
    for (unsigned int j=0; j<n_dof; ++j)
      (*scattering_src)[i](j) = 0.;

  vector<CELL*>::iterator cell(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
  for (; cell<cell_end; ++cell)
  {
    const unsigned int j_min((*cell)->Get_first_dof());
    const unsigned int j_max((*cell)->Get_last_dof());
    FINITE_ELEMENT const* const fe((*cell)->Get_fe());
    unsigned int dof_per_cell(fe->Get_dof_per_cell());
    Teuchos::SerialDenseVector<int,double> x_cell(dof_per_cell);
    Teuchos::SerialDenseVector<int,double>* scat_src_cell(
        teuchos_vector[dof_per_cell]);
    Teuchos::SerialDenseMatrix<int,double> const* const mass_matrix(
        fe->Get_mass_matrix());
    Teuchos::BLAS<int,double> blas;
    for (unsigned int i=0; i<n_mom; ++i)
    {
      for (unsigned int j=j_min; j<j_max; ++j)
        x_cell.values()[j-j_min] = x[0][i*n_dof+j];
      blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,
          (*cell)->Get_sigma_s(group,lvl,quad->Get_l(i)),mass_matrix->values(),
          mass_matrix->stride(),x_cell.values(),1,0.,scat_src_cell->values(),1);
      for (unsigned int j=j_min; j<j_max; ++j)
        (*scattering_src)[i].values()[j] += scat_src_cell->values()[j-j_min];
    }
  }
}

void TRANSPORT_OPERATOR::Sweep(Epetra_MultiVector &flux_moments,bool rhs) const
{
  const unsigned int n_cells(dof_handler->Get_n_cells());
  const unsigned int n_dir(quad->Get_n_dir());
  const unsigned int n_mom(quad->Get_n_mom());
  Teuchos::SerialDenseMatrix<int,double> const* const M2D(quad->Get_M2D());
  Teuchos::SerialDenseMatrix<int,double> const* const D2M(quad->Get_D2M());
  Epetra_Map psi_map(n_dof,0,*comm);
  Teuchos::BLAS<int,double> blas;
  Teuchos::LAPACK<int,double> lapack;
  // Clear flux_moments
  for (unsigned int i=0; i<n_dof*n_mom; ++i)
    flux_moments[0][i] = 0.;
  // Loop on the direction
  for (unsigned int idir=0; idir<n_dir; ++idir)
  { 
    Epetra_MultiVector psi(psi_map,1);
    // Get the direction
    Teuchos::SerialDenseVector<int,double> omega(quad->Get_omega_2d(idir));
    const double omega_0(omega(0));
    const double omega_1(omega(1));
    
    // Sweep on the spatial cells
    ui_vector const* const sweep_order(dof_handler->Get_sweep_order(lvl,idir));
    for (unsigned int i=0; i<n_cells; ++i)
    {
      CELL* cell(dof_handler->Get_cell((*sweep_order)[i]));
      FINITE_ELEMENT const* const fe(cell->Get_fe());
      const unsigned int dof_per_cell(fe->Get_dof_per_cell());
      const unsigned int dof_per_cell_square(dof_per_cell*dof_per_cell);
      const unsigned int offset(cell->Get_first_dof());
      Teuchos::SerialDenseVector<int,double> b(dof_per_cell);
      Teuchos::SerialDenseMatrix<int,double> const* const mass_matrix(
          fe->Get_mass_matrix());
      Teuchos::SerialDenseMatrix<int,double> const* const x_grad_matrix(
          fe->Get_x_grad_matrix());
      Teuchos::SerialDenseMatrix<int,double> const* const y_grad_matrix(
          fe->Get_y_grad_matrix());
      Teuchos::SerialDenseMatrix<int,double> A(*(fe->Get_mass_matrix()));
      
      // Volumetric term of the lhs : 
      // -omega_x * x_grad_matrix - omega_y *y_grad_matrix + sigma_t mass
      A.scale(cell->Get_sigma_t(group,lvl));
      for (unsigned int i=0; i<dof_per_cell_square; ++i)
        A.values()[i] += -omega_0*x_grad_matrix->values()[i]-omega_1*
          y_grad_matrix->values()[i];

      // Volumetric term of the rhs
      for (unsigned int mom=0; mom<n_mom; ++mom)
      {
        const double m2d((*M2D)(idir,mom));
        double* scat_src_val((*scattering_src)[mom].values());
        for (unsigned int j=0; j<dof_per_cell; ++j)
          b.values()[j] += m2d*scat_src_val[offset+j];
      }
      if (rhs==true)
      {
        // Divide the source by the sum of the weights so the input source is 
        // easier to set
        const double src(cell->Get_source(group)/param->Get_weight_sum());
        for (unsigned int k=0; k<dof_per_cell; ++k)
        {
          double const* mass_matrix_k((*mass_matrix)[k]);
          for (unsigned int j=0; j<dof_per_cell; ++j)
              b.values()[j] += src*mass_matrix_k[j];
        }
      }
      // Surfacic terms
      bool reflective_b(false);
      unsigned int edge_lid(0);
      vector<EDGE*>::iterator cell_edge(cell->Get_cell_edges_begin());
      vector<EDGE*>::iterator cell_edge_end(cell->Get_cell_edges_end());
      for (; cell_edge<cell_edge_end; ++cell_edge)
      {
        unsigned int index_cell(0);
        if ((*cell_edge)->Get_cell_index(0)!=cell->Get_id())
          index_cell = 1;
        Teuchos::SerialDenseVector<int,double> const* const external_normal(
            (*cell_edge)->Get_external_normal(index_cell));
        const double n_dot_omega(omega.dot(*external_normal));
        if ((*cell_edge)->Is_interior()==false)
          if ((*cell_edge)->Get_bc_type()==reflective)
            reflective_b = true;
        if (n_dot_omega<0.)
        {
          // Upwind
          // Check if the edge is on the border and if rhs is true
          if ((*cell_edge)->Is_interior()==true)
          {
            CELL* upwind_cell(NULL);
            if ((*cell_edge)->Get_cell_index(0)==cell->Get_id())
              upwind_cell = dof_handler->Get_cell((*cell_edge)->Get_cell_index(1));
            else
              upwind_cell = dof_handler->Get_cell((*cell_edge)->Get_cell_index(0));
            FINITE_ELEMENT const* const upwind_fe(upwind_cell->Get_fe());
            const unsigned int j_min(upwind_cell->Get_first_dof());
            const unsigned int j_max(upwind_cell->Get_last_dof());
            Teuchos::SerialDenseVector<int,double>* psi_cell(
                teuchos_vector[upwind_fe->Get_dof_per_cell()]);
            Teuchos::SerialDenseMatrix<int,double> const* const upwind(
                fe->Get_upwind_matrix(edge_lid));
            for (unsigned int j=j_min; j<j_max; ++j)
              psi_cell->values()[j-j_min] = psi[0][j];

            blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,j_max-j_min,-n_dot_omega,
                upwind->values(),upwind->stride(),psi_cell->values(),1,1.,
                b.values(),1);
          }                                                                         
          else
          {
            Teuchos::SerialDenseMatrix<int,double> const* const downwind(
                fe->Get_downwind_matrix(edge_lid));
            if ((rhs==true) && ((*cell_edge)->Get_bc_type()!=reflective))
            {
              double inc_flux_norm(0.);
              if ((*cell_edge)->Get_edge_type()==bottom_boundary)
              {
                if ((((*cell_edge)->Get_bc_type()==most_normal) &&
                      (dof_handler->Is_most_normal_bottom(lvl,idir))) ||
                    ((*cell_edge)->Get_bc_type()==isotropic))
                  inc_flux_norm = param->Get_inc_bottom();
              }
              if ((*cell_edge)->Get_edge_type()==right_boundary) 
              {
                if ((((*cell_edge)->Get_bc_type()==most_normal) &&
                      (dof_handler->Is_most_normal_right(lvl,idir))) ||
                    ((*cell_edge)->Get_bc_type()==isotropic))
                  inc_flux_norm = param->Get_inc_right();
              }
              if ((*cell_edge)->Get_edge_type()==top_boundary) 
              {
                if ((((*cell_edge)->Get_bc_type()==most_normal) &&
                      (dof_handler->Is_most_normal_top(lvl,idir))) ||
                    ((*cell_edge)->Get_bc_type()==isotropic))
                  inc_flux_norm = param->Get_inc_top();
              }
              if ((*cell_edge)->Get_edge_type()==left_boundary) 
              {
                if ((((*cell_edge)->Get_bc_type()==most_normal) &&
                      (dof_handler->Is_most_normal_left(lvl,idir))) ||
                    ((*cell_edge)->Get_bc_type()==isotropic))
                  inc_flux_norm = param->Get_inc_left();
              }
              inc_flux_norm /= param->Get_weight_sum();
              Teuchos::SerialDenseVector<int,double> inc_flux(dof_per_cell);
              inc_flux.putScalar(inc_flux_norm);
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,-n_dot_omega,
                  downwind->values(),downwind->stride(),inc_flux.values(),1,1.,
                  b.values(),1);
            }
            if ((*cell_edge)->Get_bc_type()==reflective)
            {
              Teuchos::SerialDenseVector<int,double> inc_flux(
                  Get_saf(idir,n_dir,n_mom,dof_per_cell,flux_moments,cell,
                    *cell_edge));
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,-n_dot_omega,
                  downwind->values(),downwind->stride(),inc_flux.values(),
                  1,1.,b.values(),1);
            }
          }
        }
        else
        {
          // Downwind
          Teuchos::SerialDenseMatrix<int,double> const* const downwind(
              fe->Get_downwind_matrix(edge_lid));
          for (unsigned int i=0; i<dof_per_cell_square; ++i)
            A.values()[i] += n_dot_omega*downwind->values()[i];
        }
        ++edge_lid;
      }
      int info;
      int ipiv[dof_per_cell];
      // Compute an LU factorization of a general M-by-N matrix using partial
      // pivoting with row interchanges
      lapack.GETRF(dof_per_cell,dof_per_cell,A.values(),A.stride(),ipiv,&info);
      // Solve a system of linear equations Ax=b with a general N-by-N matrix A
      // usinf the LU factorization computed by GETRF
      lapack.GETRS('N',dof_per_cell,1,A.values(),A.stride(),ipiv,
          b.values(),b.stride(),&info);

      const unsigned int j_min(cell->Get_first_dof());
      const unsigned int j_max(cell->Get_last_dof());
      for (unsigned int j=j_min; j<j_max; ++j)
        psi[0][j] = b.values()[j-j_min];
      
      if (reflective_b==true)
        Store_saf(psi,flux_moments,cell,idir,n_mom,dof_per_cell);
    }
    // Update scalar flux
    for (unsigned int mom=0; mom<n_mom; ++mom)
    {
      const unsigned int offset(mom*n_dof);
      const double d2m((*D2M)(mom,idir));
      for (unsigned int j=0; j<n_dof; ++j)
        flux_moments[0][j+offset] += d2m*psi[0][j];
    }
  }
  // Update the significant angular fluxes
  if (dof_handler->Get_n_sf_per_dir()>0)
    Update_saf(flux_moments,n_dir);
}

Teuchos::SerialDenseVector<int,double> TRANSPORT_OPERATOR::Get_saf(
    unsigned int idir,unsigned int n_dir,unsigned int n_mom,
    unsigned int dof_per_cell,Epetra_MultiVector &flux_moments,
    CELL const* const cell,EDGE const* const edge) const
{
  const unsigned int n_dir_quad(n_dir/4);
  unsigned int reflec_dir(idir);
  unsigned int offset(0);
  Teuchos::SerialDenseVector<int,double> values(dof_per_cell);

  if (edge->Get_edge_type()==left_boundary)
  {
    if (idir<n_dir_quad)
      reflec_dir = idir+3*n_dir_quad;
    else
    {
      assert(idir<2*n_dir_quad);
      reflec_dir = idir+n_dir_quad;
    }
  }
  if (edge->Get_edge_type()==right_boundary)
  {
    if ((idir>=2*n_dir_quad) && (idir<3.*n_dir_quad))
      reflec_dir = idir-n_dir_quad;
    else
    {
      assert(idir>=3*n_dir_quad);
      reflec_dir = idir-3*n_dir_quad;
    }
  }
  if (edge->Get_edge_type()==top_boundary)
  {
    if ((idir>=n_dir_quad) && (idir<2*n_dir_quad))
      reflec_dir = idir-n_dir_quad;
    else
    {
      assert((idir>=2*n_dir_quad) && (idir<3*n_dir_quad));
      reflec_dir = idir+n_dir_quad;
    }
  }
  if (edge->Get_edge_type()==bottom_boundary)
  {
    if (idir<n_dir_quad)
      reflec_dir = idir+n_dir_quad;
    else
    {
      assert(idir>=3*n_dir_quad);
      reflec_dir = idir-n_dir_quad;
    }
  }

  offset = dof_handler->Get_n_dof()*n_mom+
    dof_handler->Get_saf_map_reflective_dof(cell->Get_id())+
    reflec_dir*dof_handler->Get_n_sf_per_dir();
  
  for (unsigned int i=0; i<dof_per_cell; ++i)
    values.values()[i] = flux_moments[0][offset+i];

  return values;
}

void TRANSPORT_OPERATOR::Store_saf(Epetra_MultiVector const &psi,
    Epetra_MultiVector &flux_moments,CELL const* const cell,unsigned int idir,
    unsigned int n_mom,unsigned int dof_per_cell) const
{
  unsigned int offset_1(dof_handler->Get_saf_map_reflective_dof(cell->Get_id())+
      idir*dof_handler->Get_n_sf_per_dir());
  unsigned int offset_2(dof_handler->Get_saf_map_dof(cell->Get_id()));

  for (unsigned int i=0; i<dof_per_cell; ++i)
    (*saf)[0][offset_1+i] = psi[0][offset_2+i];
}

void TRANSPORT_OPERATOR::Update_saf(Epetra_MultiVector &flux_moments,
    unsigned int n_dir) const
{
  const unsigned int offset(dof_handler->Get_n_dof()*quad->Get_n_mom());
  const unsigned int i_max(n_dir*dof_handler->Get_n_sf_per_dir());
  for (unsigned int i=0; i<i_max; ++i)
    flux_moments[0][offset+i] = (*saf)[0][i];
}

char const* TRANSPORT_OPERATOR::Label() const
{
  return "transport_operator";
}

bool TRANSPORT_OPERATOR::UseTranspose() const
{
  return false;
}

bool TRANSPORT_OPERATOR::HasNormInf() const
{
  return false;
}

Epetra_Comm const& TRANSPORT_OPERATOR::Comm() const
{
  return *comm;
}

Epetra_Map const& TRANSPORT_OPERATOR::OperatorDomainMap() const
{
  return *flux_moments_map;
}

Epetra_Map const& TRANSPORT_OPERATOR::OperatorRangeMap() const
{
  return *flux_moments_map;
}

Epetra_MultiVector TRANSPORT_OPERATOR::Restrict_vector(Epetra_MultiVector &x) const
{
  const unsigned int i_max(dof_handler->Get_n_dof()*quad->Get_n_mom());
  Epetra_MultiVector restriction(*flux_moments_map,1);
  if (param->Get_galerkin()==true)
  {
    // With Galerkin some moments needs to be skipped because of the selection
    // rules
    double sn(-1.+sqrt(1.+2.*quad->Get_n_dir()));
    const unsigned int copy((pow(sn,2.)+sn)/2.*dof_handler->Get_n_dof());
    const unsigned int skip((sn/2.+1.)*dof_handler->Get_n_dof());
    for (unsigned int i=0; i<copy; ++i)
      restriction[0][i] = x[0][i];
    for (unsigned int i=copy; i<i_max; ++i)
      restriction[0][i] = x[0][skip+i];
  }
  else
  {
    for (unsigned int i=0; i<i_max; ++i)
      restriction[0][i] = x[0][i];
  }
   
  return restriction;
}

void TRANSPORT_OPERATOR::Project_vector(Epetra_MultiVector &x,Epetra_MultiVector &y)
  const
{
  const unsigned int y_size(y.MyLength());
  if (param->Get_galerkin()==true)
  {
    // With Galerkin some moments needs to be skipped because of the selection
    // rules
    double n_dir(y_size/dof_handler->Get_n_dof());
    double sn(-1.+sqrt(1.+2.*n_dir));
    const unsigned int common((pow(sn,2.)+sn)/2.*dof_handler->Get_n_dof());
    const unsigned int skip((sn/2.+1.)*dof_handler->Get_n_dof());
    for (unsigned int i=0; i<common; ++i)
      x[0][i] += y[0][i];
    for (unsigned int i=common; i<y_size; ++i)
      x[0][skip+i] += y[0][i];
  }
  else
  {
    for (unsigned int i=0; i<y_size; ++i)
      x[0][i] += y[0][i];
  }
}

void TRANSPORT_OPERATOR::Set_group(unsigned int g)
{
  group = g;
  if (param->Get_mip()==true)
    precond->Set_group(g);
}
