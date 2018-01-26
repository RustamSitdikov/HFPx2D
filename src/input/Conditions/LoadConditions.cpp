//
// Created by Federico Ciardo on 13/11/17.
//

#include "LoadConditions.h"

namespace hfp2d {

InSituStress loadConditions(
    Mesh &theLoadedMesh, const il::String &inputFileName,
    const il::MapArray<il::String, il::Dynamic> &ConditionsMap) {

  double far_field_vert_stress =
      findDouble("far_field_vertical_stress", ConditionsMap, inputFileName);
  double far_field_horiz_stress =
      findDouble("far_field_horizontal_stress", ConditionsMap, inputFileName);
  double amb_pore_pressure =
      findDouble("ambient_pore_pressure", ConditionsMap, inputFileName);

  // creation of the array/vectors
  il::Array2D<double> insitu_stress_distribution(
      2 * theLoadedMesh.numberOfElts(), 2);
  il::Array<double> amb_pore_pressure_distribution(
          theLoadedMesh.numberPressDofs());

  il::Array2D<il::int_t> dof_single_dd{
          theLoadedMesh.numberOfElts(),
      (theLoadedMesh.interpolationOrder() + 1), 0};
  for (il::int_t i = 0; i < theLoadedMesh.numberOfElts(); i++) {
    for (il::int_t j = 0; j < 1 * (theLoadedMesh.interpolationOrder() + 1);
         j++) {
      dof_single_dd(i, j) = i * 1 * (theLoadedMesh.interpolationOrder() + 1) + j;
    }
  }

  for (int elmt_k = 0; elmt_k < theLoadedMesh.numberOfElts(); ++elmt_k) {
    hfp2d::SegmentData mysege = theLoadedMesh.getElementData(elmt_k);

    for (int coll_k = 0; coll_k < theLoadedMesh.interpolationOrder() + 1;
         ++coll_k) {

      insitu_stress_distribution(dof_single_dd(elmt_k, coll_k), 0) =
          (far_field_vert_stress * sin(mysege.theta())) +
          (far_field_horiz_stress * cos(mysege.theta()));
      insitu_stress_distribution(dof_single_dd(elmt_k, coll_k), 1) =
          (far_field_vert_stress * cos(mysege.theta())) +
          (far_field_horiz_stress * sin(mysege.theta()));

    }
  }

  // Simple and homogeneous pore pressure distribution at nodal points within
  // the crack at ambient conditions
  for (il::int_t k = 0; k < amb_pore_pressure_distribution.size(); ++k) {
    amb_pore_pressure_distribution[k] = amb_pore_pressure;
  }

  // Call the constructor for Loading Conditions
  InSituStress BackgroundLoadingConditions(insitu_stress_distribution,
                                 amb_pore_pressure_distribution);

  return BackgroundLoadingConditions;
};
}