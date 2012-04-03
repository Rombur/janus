#ifndef _TRANSPORT_OPERATOR_HH_
#define _TRANSPORT_OPERATOR_HH_

#include <cassert>
#include "gsl_math.h"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_BLAS_types.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "DOF_HANDLER.hh"
#include "MIP.hh"
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
 *
 */
class TRANSPORT_OPERATOR : public Epetra_Operator
{
  public :
    TRANSPORT_OPERATOR(DOF_HANDLER* dof,PARAMETERS const* param,
        QUADRATURE const* quad,Epetra_Comm const* comm,
        Epetra_Map const* flux_moments_map);

    ~TRANSPORT_OPERATOR();

    /// Return the result of the transport operator applied to x in y. Return
    /// 0 if successful.
    int Apply(Epetra_MultiVector const &x,Epetra_MultiVector &y) const;

    /// Compute the scattering given a flux.
    void Compute_scattering_source(Epetra_MultiVector const &x) const;

    /// Perform a sweep. If rhs is false, the surfacic and the volumetric
    /// sources are not included in the sweep.
    /// @todo The sweep can be optimized to use less memory (only keep the
    /// front wave).
    void Sweep(Epetra_MultiVector &flx_moments,bool rhs=false) const;

    /// This method is not implemented.
    int SetUseTranspose(bool UseTranspose) {};

    /// This method is not implemented.
    int ApplyInverse(Epetra_MultiVector const &X, Epetra_MultiVector &Y) const {};

    /// This method is not implemented.
    double NormInf() const {};

    /// Return a character string describing the operator.
    char const* Label() const;

    /// Return the UseTranspose setting (always false).
    bool UseTranspose() const;

    /// Return true if the this object can provide an approximate Inf-norm,
    /// false otherwise.
    bool HasNormInf() const;

    /// Return a pointer to the Epetra_Comm communicator associated with this
    /// operator.
    Epetra_Comm const& Comm() const;

    /// Return the Epetra_Map object associated with the domain of this
    /// operator.
    Epetra_Map const& OperatorDomainMap() const;

    /// Return the Epetra_Map object associated with the range of this
    /// operator.
    Epetra_Map const& OperatorRangeMap() const;

  private :
    /// Return a Teuchos vector with the value of the significant flux for
    /// a given cell and a given direction.
    Teuchos::SerialDenseVector<int,double> Get_saf(unsigned int idir,
        unsigned int n_dir,unsigned int dof_per_cell,Epetra_MultiVector &flx_moments,
        CELL const* const cell,EDGE const* const edge) const;

    /// Store the Significant Angular Fluxes of a given cell and direction.
    void Store_saf(Epetra_MultiVector const &psi,Epetra_MultiVector &flx_moments,
        CELL const* const cell,unsigned int idir,unsigned int n_dir,
        unsigned int dof_per_cell) const;

    /// Temporary variable.
    unsigned int lvl;
    /// Epetra communicator.
    Epetra_Comm const* comm;
    /// Epetra map
    Epetra_Map const* flx_moments_map;
    /// Scattering source for each moment.
    /// @todo Change the Teuchos::SerialDenseVector to Epetra_MultiVector.
    vector<Teuchos::SerialDenseVector<int,double> >* scattering_src;
    /// Pointer to the dof_handler object. (DOF_HANDLER is not const because
    /// we need to return an iterator on the cells).
    DOF_HANDLER* dof_handler;
    /// Pointer to the MIP DSA.
    MIP* precond;
    /// Pointer to the parameters object.
    PARAMETERS const* param;
    /// Pointer to the quadrature object.
    QUADRATURE const* quad;
};

#endif
