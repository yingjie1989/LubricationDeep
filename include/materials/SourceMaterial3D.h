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

#ifndef SOURCEMATERIAL3D_H
#define SOURCEMATERIAL3D_H

#include "Material.h"


//Forward Declarations
class SourceMaterial3D;

template<>
InputParameters validParams<SourceMaterial3D>();

/**
 * Example material class that defines a few properties.
 */
class SourceMaterial3D : public Material
{
public:
  SourceMaterial3D(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

private:
  /**
   * This is the member reference that will hold the computed values
   * for the Real value property in this class.
   */
  MaterialProperty<Real> & _diffx;
  MaterialProperty<Real> & _diffy;

 
  Real _coeff;
  Real _dt;
  Real _mu;
  Real _component;
  /**
   * Computed values for the Gradient value property in this class.
   */
  //MaterialProperty<RealGradient> & _convection_velocity;

  /**
   * This is the member reference that will hold the gradient
   * of the coupled variable
   */
 // const VariableGradient & _diffusion_gradient;

  /**
   * This object returns a piecewise linear function based an a series
   * of points and their corresponding values
   */
  //LinearInterpolation _piecewise_func;
};

#endif //SOURCEMATERIAL3D_H
