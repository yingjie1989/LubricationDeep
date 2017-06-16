/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTRACTIONFORMNODIMEN_H
#define INSMOMENTUMTRACTIONFORMNODIMEN_H

#include "INSMomentumBaseALE.h"

// Forward Declarations
class INSMomentumTractionFormNoDimen;

template<>
InputParameters validParams<INSMomentumTractionFormNoDimen>();

/**
 * This class computes momentum equation residual and Jacobian viscous
 * contributions for the "traction" form of the governing equations.
 */
class INSMomentumTractionFormNoDimen : public INSMomentumBaseALE
{
public:
  INSMomentumTractionFormNoDimen(const InputParameters & parameters);

  virtual ~INSMomentumTractionFormNoDimen(){}

 

protected:
  virtual Real computeQpResidualViscousPart() override;
  virtual Real computeQpJacobianViscousPart() override;
  virtual Real computeQpOffDiagJacobianViscousPart(unsigned jvar) override;

  unsigned _dim;
  std::vector<Real> _coeff;

};


#endif
