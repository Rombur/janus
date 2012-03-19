#ifndef _TRANSPORT_OPERATOR_HH_
#define _TRANSPORT_OPERATOR_HH_

#include <cassert>
#include "Epetra_FEVector.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "DOF_HANDLER.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"

using namespace std;

typedef vector<unsigned int> ui_vector;

/**
 * This class derives from Epetra_Operator and implement the function Apply
 * needed by the AztecOO solvers and the compute the sweep. 
 * @todo The reflective boundary can be optimized. Right now it uses to much
 * memory.
 * @todo Reflective boundary with angular multigrid.
 */
class TRANSPORT_OPERATOR : public Epetra_Operator
{
  public :
    TRANSPORT_OPERATOR(DOF_HANDLER const* dof,PARAMETERS const* param,
        QUADRATURE* quad);

    ~TRANSPORT_OPERATOR();

    /// Return the result of the transport operator applied to x in y. Return
    /// 0 if successful.
    int Apply(Epetra_FEVector const &x,Epetra_FEVector &y);

    /// Compute the scattering given a flux.
    void Compute_scattering_source(Epetra_FEVector const &x);

  private :
    /// Return a Teuchos vector with the value of the significant flux for
    /// a given cell and a given direction.
    Teuchos::SerialDenseVector<int,double> Get_saf(unsigned int idir,
        unsigned int n_dir,unsigned int dof_per_cell,Epretra_FEVector &flx_moments,
        CELL const* const cell,EDGE const* const edge) const;

    /// Store the Significant Angular Fluxes of a given cell and direction.
    void Store_saf(Epetra_FEVector const &psi,CELL const* const cell,
    unsigned int idir,unsigned int n_dir,unsigned int n_dof);

    /// Temporary variable.
    unsigned int lvl;
    vector<Teuchos::SerialDenseVector<int,double> > scattering_src;
    /// Pointer to the dof_handler object.
    DOF_HANDLER const* dof_handler;
    /// Pointer to the parameters object.
    PARAMETERS const* param;
    /// Pointer to the quadrature object.
    QUADRATURE* quad;
};

#endif
