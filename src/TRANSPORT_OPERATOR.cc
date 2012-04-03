#include "TRANSPORT_OPERATOR.hh"

TRANSPORT_OPERATOR::TRANSPORT_OPERATOR(DOF_HANDLER* dof,
    PARAMETERS const* param, QUADRATURE const* quad,Epetra_Comm const* comm,
    Epetra_Map const* flux_moments_map) :
  Epetra_Operator(),
  lvl(0),
  comm(comm),
  flx_moments_map(flux_moments_map),
  scattering_src(NULL),
  dof_handler(dof),
  precond(NULL),
  param(param),
  quad(quad)
{
  scattering_src = new vector<Teuchos::SerialDenseVector<int,double> >
    (quad->Get_n_mom(),Teuchos::SerialDenseVector<int,double> 
     (dof_handler->Get_n_dof()));
  if (param->Get_mip()==true)
    precond = new MIP (dof,param,quad,comm);
}

TRANSPORT_OPERATOR::~TRANSPORT_OPERATOR()
{
  if (precond!=NULL)
  {
    delete precond;
    precond = NULL;
  }

  if (scattering_src!=NULL)
  {
    delete scattering_src;
    scattering_src = NULL;
  }
}

int TRANSPORT_OPERATOR::Apply(Epetra_MultiVector const &x,Epetra_MultiVector &y) const
{
  const unsigned int n_dof(dof_handler->Get_n_dof());
  y = x;

  if (param->Get_multigrid()==true)
  {
  }
  else
  {
    if (param->Get_mip()==true)
      precond->Solve(y);

    // Compute the scattering source
    Compute_scattering_source(y);
    Sweep(y);

    for (unsigned int i=0; i<n_dof; ++i)
      y[0][i] = x[0][i]-y[0][i];
  }  

  return 0;
}

void TRANSPORT_OPERATOR::Compute_scattering_source(Epetra_MultiVector const &x) const
{
  const unsigned int n_mom(quad->Get_n_mom());
  const unsigned int n_dof(dof_handler->Get_n_dof());
  // Reinitialize the scattering source
  for (unsigned int i=0; i<n_mom; ++i)
    for (unsigned int j=0; j<n_dof; ++j)
      (*scattering_src)[i](j) = 0.;

  vector<CELL*>::iterator cell(dof_handler->Get_mesh_begin());
  vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
  for (; cell<cell_end; ++cell)
  {
    for (unsigned int i=0; i<n_mom; ++i)
    {
      const unsigned int j_min((*cell)->Get_first_dof());
      const unsigned int j_max((*cell)->Get_last_dof());
      FINITE_ELEMENT const* const fe((*cell)->Get_fe());
      unsigned int dof_per_cell(fe->Get_dof_per_cell());
      Teuchos::SerialDenseVector<int,double> x_cell(dof_per_cell);
      Teuchos::SerialDenseVector<int,double> scat_src_cell(dof_per_cell);
      Teuchos::SerialDenseMatrix<int,double> const* const mass_matrix(
          fe->Get_mass_matrix());
      Teuchos::BLAS<int,double> blas;
      for (unsigned int j=j_min; j<j_max; ++j)
        x_cell(j-j_min) = x[0][j];
      blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,
          (*cell)->Get_sigma_s(lvl,i),mass_matrix->values(),
          mass_matrix->stride(),x_cell.values(),1,0.,scat_src_cell.values(),1);
      for (unsigned int j=j_min; j<j_max; ++j)
        (*scattering_src)[i](j) += scat_src_cell(j-j_min);
    }
  }
}

void TRANSPORT_OPERATOR::Sweep(Epetra_MultiVector &flx_moments,bool rhs) const
{
  const unsigned int n_cells(dof_handler->Get_n_cells());
  const unsigned int n_dir(quad->Get_n_dir());
  const unsigned int n_mom(quad->Get_n_mom());
  const unsigned int n_dof(dof_handler->Get_n_dof());
  Teuchos::SerialDenseMatrix<int,double> const* const M2D(quad->Get_M2D());
  Teuchos::SerialDenseMatrix<int,double> const* const D2M(quad->Get_D2M());
  Epetra_Map psi_map(n_dof,0,*comm);
  Teuchos::BLAS<int,double> blas;
  Teuchos::LAPACK<int,double> lapack;
  // Clear flx_moments
  flx_moments.Scale(0.);
  // Loop on the direction
  for (unsigned int idir=0; idir<n_dir; ++idir)
  {
    Epetra_MultiVector psi(psi_map,1);
    // Get the direction
    Teuchos::SerialDenseVector<int,double> omega(quad->Get_omega_2d(idir));
    
    // Sweep on the spatial cells
    ui_vector const* const sweep_order(dof_handler->Get_sweep_order(lvl,idir));
    for (unsigned int i=0; i<n_cells; ++i)
    {
      CELL* cell(dof_handler->Get_cell((*sweep_order)[i]));
      FINITE_ELEMENT const* const fe(cell->Get_fe());
      const unsigned int dof_per_cell(fe->Get_dof_per_cell());
      const unsigned int offset(cell->Get_first_dof());
      Teuchos::SerialDenseVector<int,double> b(dof_per_cell);
      Teuchos::SerialDenseMatrix<int,double> A(*(fe->Get_x_grad_matrix()));
      Teuchos::SerialDenseMatrix<int,double> tmp(*(fe->Get_y_grad_matrix()));
      Teuchos::SerialDenseMatrix<int,double> const* const mass_matrix(
          fe->Get_mass_matrix());
      
      // Volumetric term of the lhs : 
      // -omega_x * x_grad_matrix - omega_y *y_grad_matrix + sigma_t mass
      A *= -omega(0);
      tmp *= -omega(1);
      A += tmp;
      tmp = *mass_matrix;
      tmp *= cell->Get_sigma_t(lvl);
      A += tmp;

      // Volumetric term of the rhs
      for (unsigned int mom=0; mom<n_mom; ++mom)
        for (unsigned int j=0; j<dof_per_cell; ++j)
          b(j) += (*M2D)(idir,mom)*(*scattering_src)[mom](offset+j);
      if (rhs==true)
      {
        // Divide the source by 4 PI so the input source is easier to set
        for (unsigned int j=0; j<dof_per_cell; ++j)
          for (unsigned int k=0; k<dof_per_cell; ++k)
            b(j) += cell->Get_source()*(*mass_matrix)(j,k)/(4.*M_PI);
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
        if ((*cell_edge)->Is_reflective()==true)
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
            const unsigned int j_min(upwind_cell->Get_first_dof());
            const unsigned int j_max(upwind_cell->Get_last_dof());
            Teuchos::SerialDenseVector<int,double> psi_cell(fe->Get_dof_per_cell());
            Teuchos::SerialDenseMatrix<int,double> const* const upwind(
                fe->Get_upwind_matrix(edge_lid));
            for (unsigned int j=j_min; j<j_max; ++j)
              psi_cell(j-j_min) = psi[0][j];
            blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,-n_dot_omega,
                upwind->values(),upwind->stride(),psi_cell.values(),1,1.,
                b.values(),1);
          }                                                                         
          else
          {
            Teuchos::SerialDenseMatrix<int,double> const* const downwind(
                fe->Get_downwind_matrix(edge_lid));
            if ((rhs==true) && ((*cell_edge)->Is_reflective()==false))
            {
              double inc_flux_norm(0.);
              if (((*cell_edge)->Get_edge_type()==bottom_boundary) &&
                  (dof_handler->Is_most_normal_bottom(lvl,(*cell_edge)->Get_gid())))
                inc_flux_norm = param->Get_inc_bottom();
              if (((*cell_edge)->Get_edge_type()==right_boundary) &&
                  (dof_handler->Is_most_normal_right(lvl,(*cell_edge)->Get_gid())))
                inc_flux_norm = param->Get_inc_right();
              if (((*cell_edge)->Get_edge_type()==top_boundary) &&
                  (dof_handler->Is_most_normal_top(lvl,(*cell_edge)->Get_gid())))
                inc_flux_norm = param->Get_inc_top();
              if (((*cell_edge)->Get_edge_type()==left_boundary) &&
                  (dof_handler->Is_most_normal_left(lvl,(*cell_edge)->Get_gid())))
                inc_flux_norm = param->Get_inc_left();
              for (unsigned int j=0; j<dof_per_cell; ++j)
                for (unsigned int k=0; k<dof_per_cell; ++k)
                  b(j) -= n_dot_omega*inc_flux_norm*(*downwind)(j,k);
            }
            if ((*cell_edge)->Is_reflective()==true)
            {
              Teuchos::SerialDenseVector<int,double> inc_flux(
                  Get_saf(idir,n_dir,dof_per_cell,flx_moments,cell,*cell_edge));
              blas.GEMV(Teuchos::NO_TRANS,dof_per_cell,dof_per_cell,-1.,
                  downwind->values(),downwind->stride(),inc_flux.values(),
                  1,-1.,b.values(),1);
            }
          }
        }
        else
        {
          // Downwind
          Teuchos::SerialDenseMatrix<int,double> tmp(
              *(fe->Get_downwind_matrix(edge_lid)));
          tmp *= n_dot_omega;
          A += tmp;
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
        psi[0][j] = b(j-j_min);
      if (reflective_b==true)
        Store_saf(psi,flx_moments,cell,idir,n_dir,dof_per_cell);
    }
    // Update scalar flux
    for (unsigned int mom=0; mom<n_mom; ++mom)
    {
      const unsigned int offset(mom*n_dof);
      for (unsigned int j=0; j<n_dof; ++j)
        flx_moments[0][j+offset] += (*D2M)(mom,idir)*psi[0][j];
    }
  }
}

Teuchos::SerialDenseVector<int,double> TRANSPORT_OPERATOR::Get_saf(
    unsigned int idir,unsigned int n_dir,unsigned int dof_per_cell,
    Epetra_MultiVector &flx_moments,CELL const* const cell,
    EDGE const* const edge) const
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
    if ((idir>n_dir_quad) && (idir<2*n_dir_quad))
      reflec_dir = idir-n_dir_quad;
    else
    {
      assert((idir>=2*n_dir_quad) && (idir<3*n_dir_quad));
      reflec_dir = idir+n_dir_quad;
    }
  }
  if (edge->Get_edge_type()==bottom_boundary)
  {
    if (idir<n_dir)
      reflec_dir = idir+n_dir_quad;
    else
    {
      assert(idir>3*n_dir_quad);
      reflec_dir = idir-n_dir_quad;
    }
  }

  offset = dof_handler->Get_n_dof()*n_dir+dof_handler->Get_saf_map(cell->Get_id())+
    reflec_dir*dof_handler->Get_n_sf_per_dir();
  
  for (unsigned int i=0; i<dof_per_cell; ++i)
    values[i] = flx_moments[0][offset+i];
}

void TRANSPORT_OPERATOR::Store_saf(Epetra_MultiVector const &psi,
    Epetra_MultiVector &flx_moments,CELL const* const cell,unsigned int idir,
    unsigned int n_dir,unsigned int dof_per_cell) const
{
  unsigned int offset_1(dof_handler->Get_n_dof()*n_dir+
      dof_handler->Get_saf_map(cell->Get_id())+idir*dof_handler->Get_n_sf_per_dir());
  unsigned int offset_2(dof_handler->Get_saf_map(cell->Get_id()));

  for (unsigned int i=0; i<dof_per_cell; ++i)
    flx_moments[0][offset_1+i] = psi[0][offset_2+i];
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
  return *flx_moments_map;
}

Epetra_Map const& TRANSPORT_OPERATOR::OperatorRangeMap() const
{
  return *flx_moments_map;
}
