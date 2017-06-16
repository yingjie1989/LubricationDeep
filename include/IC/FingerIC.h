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

#ifndef FINGERIC_H
#define FINGERIC_H

#include "InitialCondition.h"

class FingerIC;

template<>
InputParameters validParams<FingerIC>();


class FingerIC : public InitialCondition
{
public:
  FingerIC(const InputParameters & parameters);

  virtual Real value(const Point & p);

private:
  Real _b;
  Real _A;
  Real _B; 
  bool _hIC;

  Real funcH(Real x);
};


#endif /*FingerIC_H */
