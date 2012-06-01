#ifndef _TRANSPORT_OPERATOR_HH_
#define _TRANSPORT_OPERATOR_HH_

#include <cassert>
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
 */
class TRANSPORT_OPERATOR : public Epetra_Operator
{
  public :
    /// Constructor used when the angular multigrid is not used.
    TRANSPORT_OPERATOR(DOF_HANDLER* dof,PARAMETERS const* param,
        QUADRATURE* quad,Epetra_Comm const* comm,
        Epetra_Map const* flux_moments_map);

    /// Constructor used for the angular multigrid.
    TRANSPORT_OPERATOR(DOF_HANDLER* dof,PARAMETERS const* param,
        vector<QUADRATURE*> const* quad_vector,Epetra_Comm const* comm, 
        Epetra_Map const* flux_moments_map,unsigned int level,unsigned int max_level,
        MIP* preconditioner=NULL);

    ~TRANSPORT_OPERATOR();

    /// Return the result of the transport operator applied to x in y. Return
    /// 0 if successful.
    int Apply(Epetra_MultiVector const &x,Epetra_MultiVector &y) const;

    /// Apply the angular multigrid preconditioner to a given vector.
    void Apply_preconditioner(Epetra_MultiVector &x);

    /// Compute the scattering given a flux.
    void Compute_scattering_source(Epetra_MultiVector const &x) const;

    /// Perform a sweep. If rhs is false, the surfacic and the volumetric
    /// sources are not included in the sweep.
    /// @todo The sweep can be optimized to use less memory (only keep the
    /// front wave).
    void Sweep(Epetra_MultiVector &flux_moments,bool rhs=false) const;

    /// This method is not implemented.
    int SetUseTranspose(bool UseTranspose) {return 0;};

    /// This method is not implemented.
    int ApplyInverse(Epetra_MultiVector const &X, Epetra_MultiVector &Y) const 
      {return 0;};

    /// This method is not implemented.
    double NormInf() const {return 0;};

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

    /// Restrict a given vector to size of this TRANSPORT_OPERATOR object.
    Epetra_MultiVector Restrict_vector(Epetra_MultiVector &x) const;

    /// Project and add a given vector y from a coarser TRANSPORT_OPERATOR 
    /// object to x. x and y are the same on output.
    void Project_vector(Epetra_MultiVector &x,Epetra_MultiVector &y) const;

    /// Return a pointer to the MIP DSA.
    MIP* Get_mip();

  private :
    /// Return a Teuchos vector with the value of the significant flux for
    /// a given cell and a given direction.
    Teuchos::SerialDenseVector<int,double> Get_saf(unsigned int idir,
        unsigned int n_dir,unsigned int n_mom,unsigned int dof_per_cell,
        Epetra_MultiVector &flux_moments,CELL const* const cell,
        EDGE const* const edge) const;

    /// Store the Significant Angular Fluxes of a given cell and direction.
    void Store_saf(Epetra_MultiVector const &psi,Epetra_MultiVector &flux_moments,
        CELL const* const cell,unsigned int idir,unsigned int n_mom,
        unsigned int dof_per_cell) const;

    /// Update the SAF stored in TRANSPORT_SOLVER::flux_moments at the end of the 
    /// sweep.
    void Update_saf(Epetra_MultiVector &flux_moments,unsigned int n_dir) const;

    /// Current level in the angular multigrid.
    unsigned int lvl;
    /// Maximum level in the angular multigrid.
    unsigned int max_lvl;
    /// Number of degrees of freedom associated to the problem.
    const unsigned int n_dof;
    /// Epetra communicator.
    Epetra_Comm const* comm;
    /// Epetra map associated to #flux_moments
    Epetra_Map const* flux_moments_map;
    /// Epetra map associated to the significant angular fluxes.
    Epetra_Map* saf_map;
    /// Store the Significant Angular Fluxes.
    Epetra_MultiVector* saf;
    /// To speed up the sweep the Teuchos Vector are constructed in the
    /// constructor and then reused in the sweep.
    vector<Teuchos::SerialDenseVector<int,double>* > teuchos_vector;
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
    QUADRATURE* quad;
    /// Pointer to the vector of pointers to the quadrature.
    vector<QUADRATURE*> const* quad_vector;
};

inline MIP* TRANSPORT_OPERATOR::Get_mip()
{
  return precond;
}

#endif
