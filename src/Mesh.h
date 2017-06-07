//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//

#ifndef HFPX2D_MESH_H
#define HFPX2D_MESH_H

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>

namespace hfp2d {

///// 1D mesh class
class Mesh { // class for 1D mesh of 1D segment elements ?

private:
  il::Array2D<double> node_;
  il::Array2D<int> connectivity_;

public:
  void set_values(il::Array2D<double>, il::Array2D<int>);

  double node(il::int_t k, il::int_t i);

  int connectivity(il::int_t k, il::int_t i);

  int nelts();
  int ncoor();
  il::Array2D<double> coor();
  il::Array2D<int> conn();
};

struct SegmentCharacteristic {
  double size;
  double theta; // angle w.r. to e_1
  il::StaticArray<double, 2>
      n; // unit normal to segment in global system of coordinates
  il::StaticArray<double, 2>
      s; // unit tangent to segment in global system of coordinates
  il::StaticArray<double, 2> Xmid; // segment mid points coordinates.
  il::Array2D<double>
      CollocationPoints; // collocation points in global system of coordinates
};

il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta);

SegmentCharacteristic get_segment_DD_characteristic(Mesh mesh, int const ne,
                                                    int const p);

struct LayerParameters1 {
  // id number to identify layer 1
  il::int_t id_layer1;
  // To set the layer size in the mesh
  // First element of layer 1
  il::int_t First_elem_layer1;
  // Last element of layer 1
  il::int_t Last_elem_layer1;
  // Cohesion of layer 1
  double Cohesion_layer1;
  // Peak friction coefficient of layer 1
  double Peak_fric_coeff_layer1;
  // Residual friction coefficient of layer 1
  double Resid_fric_coeff_layer1;
  // slip dw for scaling of layer 1 (see linear law in the report)
  double d_wf_layer1;
};

struct LayerParameters2 {
  // id number to identify layer 2
  il::int_t id_layer2;
  // To set the layer size in the mesh
  // First element of layer 2
  il::int_t First_elem_layer2;
  // Last element of layer 2
  il::int_t Last_elem_layer2;
  // Cohesion of layer 2
  double Cohesion_layer2;
  // Peak friction coefficient of layer 2
  double Peak_fric_coeff_layer2;
  // Residual friction coefficient
  double Resid_fric_coeff_layer2;
  // slip dw for scaling (see linear law in the report)
  double d_wf_layer2;
};

struct LayerParameters3 {
  // id number to identify layer 3
  il::int_t id_layer3;
  // To set the layer size in the mesh
  // First element of layer 3
  il::int_t First_elem_layer3;
  // Last element of layer 3
  il::int_t Last_elem_layer3;
  // Cohesion of layer 3
  double Cohesion_layer3;
  // Peak friction coefficient of layer 3
  double Peak_fric_coeff_layer3;
  // Residual friction coefficient of layer 3
  double Resid_fric_coeff_layer3;
  // slip dw for scaling of layer 3 (see linear law in the report)
  double d_wf_layer3;
};

il::Array<il::int_t> id_mesh_layers(Mesh mesh,
                                    LayerParameters1 layer_parameter1,
                                    LayerParameters2 layer_parameter2,
                                    LayerParameters3 layer_parameter3);
}

#endif // HFPX2D_MESH_H
