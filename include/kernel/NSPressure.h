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

#ifndef NSPressure_H
#define NSPressure_H

#include "Kernel.h"

// Forward Declaration
class NSPressure;

template<>
InputParameters validParams<NSPressure>();

class NSPressure : public Kernel
{
public:
  NSPressure(const InputParameters & parameters);

  const MaterialProperty<Real> & _grav;
  int _dim;
  Real _coeff;
  const VariableValue & _coupled;
  const VariableGradient & _coupled_grad;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);

};
#endif //NSPressure_H
