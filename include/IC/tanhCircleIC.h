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

#ifndef TANHCIRCLEIC_H
#define TANHCIRCLEIC_H

#include "InitialCondition.h"

class tanhCircleIC;

template<>
InputParameters validParams<tanhCircleIC>();


class tanhCircleIC : public InitialCondition
{
public:
  tanhCircleIC(const InputParameters & parameters);

  virtual Real value(const Point & p);

private:
  int _dim;
  std::vector<Real> _pos;
  Real _invalue;
  Real _outvalue;
  Real _rad;
  Real _m;
  Real _A; 
};


#endif /*tanhCircleIC_H */
