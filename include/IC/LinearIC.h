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

#ifndef LINEARIC_H
#define LINEARIC_H


#include "InitialCondition.h"

class LinearIC;

template<>
InputParameters validParams<LinearIC>();


class LinearIC : public InitialCondition
{
public:

   LinearIC(const InputParameters & parameters);

   virtual Real value(const Point & p) override;

private:
   int _dim;
   Real _coefficient;
};

#endif //LINEARIC_H
