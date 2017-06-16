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

#ifndef CAPILARYCIRCLEMATERIAL_H
#define CAPILARYCIRCLEMATERIAL_H

#include "Material.h"


//Forward Declarations
class CapilaryCircleMaterial;

template<>
InputParameters validParams<CapilaryCircleMaterial>();

/**
 * Example material class that defines a few properties.
 */
class CapilaryCircleMaterial : public Material
{
public:
  CapilaryCircleMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  /**
   * This is the member reference that will hold the computed values
   * for the Real value property in this class.
   */
  MaterialProperty<Real> & _diffusivity;
 
  Real _kappa0;
  Real _eta0;

   std::vector<Real> _coeff;
   std::vector<Real> _wavex;
   std::vector<Real> _wavey;
   Real _B;

};

#endif //CAPILARYMATERIAL_H
