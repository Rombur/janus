#ifndef _LS_HH_
#define _LS_HH_

#include <cassert>
#include "QUADRATURE.hh"

using namespace std;

/**
 * This class build the Level Symmetric quadrature.
 */

class LS : public QUADRATURE
{
  public :
    LS(unsigned int sn,unsigned int L_max,bool galerkin);

  private :
    /// Compute omega in one octant.
    void Build_octant();

    /// Compute the different omega given a set of directions.
    void Compute_omega(d_vector const &direction);
};

#endif
