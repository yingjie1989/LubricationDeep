/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTESURFACESTRESS_H
#define COMPUTESURFACESTRESS_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "RotationTensor.h"
#include "DerivativeMaterialInterface.h"
#include "Assembly.h"

/**
 * ComputeEnrichStrain is the base class for strain tensors
 */
class ComputeSurfaceStress : public DerivativeMaterialInterface<Material>
{
public:
  ComputeSurfaceStress(const InputParameters & parameters);
  virtual ~ComputeSurfaceStress() {}

protected:
  virtual void initQpStatefulProperties();

  virtual void computeQpProperties();

   MaterialProperty<RealVectorValue> & _surface_stress;
  
  /*std::vector<Real> _enrich_disp;
  std::vector<RealVectorValue> _grad_enrich_disp;

  std::vector<std::vector<MooseVariable *> > _enrich_variable;
  

  /// the current shape functions
  const VariablePhiValue & _phi;

  /// gradient of the shape function
  const VariablePhiGradient & _grad_phi;*/

  /// Coupled displacement variables
  //unsigned int _ndisp;
  //std::vector<const VariableValue *> _disp;
  //std::vector<const VariableGradient *> _grad_disp;

  //std::string _base_name;

//  MaterialProperty<RankTwoTensor> & _mechanical_strain;

//  MaterialProperty<RankTwoTensor> & _total_strain;

//  const MaterialProperty<RankTwoTensor> & _eigenstrain;

    Real _hele;
    std::vector<Real> _Re;
    Real _tol;
    unsigned _dim;
      const VariableValue & _some_var_val;
      const VariableGradient & _c_grad;
      const VariableGradient & _grad_h;
//	unsigned _h_var_number;
      const VariableSecond & _h_second;
     MooseVariable* _h_var_number;

};

#endif //COMPUTESURFACESTRESS_H
