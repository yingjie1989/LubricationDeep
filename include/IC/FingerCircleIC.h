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

#ifndef FINGERCIRCLEIC_H
#define FINGERCIRCLEIC_H

#include "InitialCondition.h"

class FingerCircleIC;

template<>
InputParameters validParams<FingerCircleIC>();


class FingerCircleIC : public InitialCondition
{
public:
  FingerCircleIC(const InputParameters & parameters);

  virtual Real value(const Point & p);

private:
  Real _b;
  Real _A;
  Real _B;
  Real _C; 
  std::vector<Real> _coeff;
  bool _hIC;

  Real funcH(Real x);
};


#endif /*FingerCircleIC_H */
