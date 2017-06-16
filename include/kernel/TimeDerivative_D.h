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

#ifndef TIMEDERIVATIVE_D_H
#define TIMEDERIVATIVE_D_H

#include "TimeKernel.h"

// Forward Declaration
class TimeDerivative_D;

template<>
InputParameters validParams<TimeDerivative_D>();

class TimeDerivative_D : public TimeKernel
{
public:
  TimeDerivative_D(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  Real _coeff;
  bool _lumping;
};

#endif //TIMEDERIVATIVE_H
