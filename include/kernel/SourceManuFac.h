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

#ifndef SOURCEMANUFAC_H
#define SOURCEMANUFAC_H

#include "Kernel.h"

class SourceManuFac;

template<>
InputParameters validParams<SourceManuFac>();


class SourceManuFac : public Kernel
{
public:
  SourceManuFac(const InputParameters & parameters);
  
  const MaterialProperty<Real> & _fmat;
 


protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int);
};


#endif /* SOURCEMANUFAC_H */
