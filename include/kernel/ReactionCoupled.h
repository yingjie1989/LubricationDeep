/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#ifndef REACTIONCOUPLED_H
#define REACTIONCOUPLED_H

#include "Kernel.h"

// Forward Declaration
class ReactionCoupled;

template<>
InputParameters validParams<ReactionCoupled>();

class ReactionCoupled : public Kernel
{
public:
  ReactionCoupled(const InputParameters & parameters);

  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _grav;
  Real _coeff;
  const VariableValue & _coupled_h;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);

};
#endif //REACTIONCOUPLED_H
