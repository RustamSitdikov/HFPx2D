//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 11.01.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <cmath>

// Inclusion from the project
#include "Friction.h"

namespace hfp2d {

// Function that returns an array that contain friction coefficient (according
// to EXPONENTIAL friction weakening law)
il::Array<double> exp_friction(LayerParameters1 &layer_parameters1,
                               LayerParameters2 &layer_parameters2,
                               LayerParameters3 &layer_parameters3,
                               const il::Array<il::int_t> &id_layers,
                               il::Array2D<int> Dofw,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - layer_parameters1 -> structure that contains all the friction parameters
  //  we need for layer 1
  //  - layer_parameters2 -> structure that contains all the friction parameters
  //  we need for layer 2
  //  - layer_parameters3 -> structure that contains all the friction parameters
  //  we need for layer 3
  //  - id_layers -> vector that contains id number of each layer
  //  - Dofw -> dof handle for piecewise linear shear or opening DDs
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> f{d.size(), 0};

  for (il::int_t i = 0; i < id_layers.size(); ++i) {

    if (id_layers[i] == layer_parameters1.id_layer1) {

      f[Dofw(i, 0)] =
          layer_parameters1.Peak_fric_coeff_layer1 -
          ((layer_parameters1.Peak_fric_coeff_layer1 -
            layer_parameters1.Resid_fric_coeff_layer1) *
           (1 - exp(-d[Dofw(i, 0)] / layer_parameters1.d_wf_layer1)));

      f[Dofw(i, 1)] =
          layer_parameters1.Peak_fric_coeff_layer1 -
          ((layer_parameters1.Peak_fric_coeff_layer1 -
            layer_parameters1.Resid_fric_coeff_layer1) *
           (1 - exp(-d[Dofw(i, 1)] / layer_parameters1.d_wf_layer1)));

    } else if (id_layers[i] == layer_parameters2.id_layer2) {

      f[Dofw(i, 0)] =
          layer_parameters2.Peak_fric_coeff_layer2 -
          ((layer_parameters2.Peak_fric_coeff_layer2 -
            layer_parameters2.Resid_fric_coeff_layer2) *
           (1 - exp(-d[Dofw(i, 0)] / layer_parameters2.d_wf_layer2)));

      f[Dofw(i, 1)] =
          layer_parameters2.Peak_fric_coeff_layer2 -
          ((layer_parameters2.Peak_fric_coeff_layer2 -
            layer_parameters2.Resid_fric_coeff_layer2) *
           (1 - exp(-d[Dofw(i, 1)] / layer_parameters2.d_wf_layer2)));

    } else if (id_layers[i] == layer_parameters3.id_layer3) {

      f[Dofw(i, 0)] =
          layer_parameters3.Peak_fric_coeff_layer3 -
          ((layer_parameters3.Peak_fric_coeff_layer3 -
            layer_parameters3.Resid_fric_coeff_layer3) *
           (1 - exp(-d[Dofw(i, 0)] / layer_parameters3.d_wf_layer3)));

      f[Dofw(i, 1)] =
          layer_parameters3.Peak_fric_coeff_layer3 -
          ((layer_parameters3.Peak_fric_coeff_layer3 -
            layer_parameters3.Resid_fric_coeff_layer3) *
           (1 - exp(-d[Dofw(i, 1)] / layer_parameters3.d_wf_layer3)));
    }
  }

  return f;
};

// Function that returns an array that contain friction coefficient (according
// to LINEAR friction weakening law)
il::Array<double> lin_friction(LayerParameters1 &layer_parameters1,
                               LayerParameters2 &layer_parameters2,
                               LayerParameters3 &layer_parameters3,
                               const il::Array<il::int_t> &id_layers,
                               il::Array2D<int> Dofw,
                               const il::Array<double> &d, il::io_t) {

  // Inputs:
  //  - layer_parameters1 -> structure that contains all the friction parameters
  //  we need for layer 1
  //  - layer_parameters2 -> structure that contains all the friction parameters
  //  we need for layer 2
  //  - layer_parameters3 -> structure that contains all the friction parameters
  //  we need for layer 3
  //  - id_layers -> vector that contains id number of each layer
  //  - Dofw -> dof handle for piecewise linear shear or opening DDs
  //  - d -> vector that contains the slip
  //  - io_t -> everything on the left of il::io_t is read-only and is not
  //    going to be mutated

  il::Array<double> f{d.size(), 0};
  double_t sl1;
  sl1 =
      layer_parameters1.Peak_fric_coeff_layer1 / layer_parameters1.d_wf_layer1;
  double_t sl2;
  sl2 =
      layer_parameters2.Peak_fric_coeff_layer2 / layer_parameters2.d_wf_layer2;
  double_t sl3;
  sl3 =
      layer_parameters3.Peak_fric_coeff_layer3 / layer_parameters3.d_wf_layer3;

  for (il::int_t i = 0; i < id_layers.size(); ++i) {

    if (id_layers[i] == layer_parameters1.id_layer1) {

      if (d[Dofw(i, 0)] < ((layer_parameters1.Peak_fric_coeff_layer1 -
                            layer_parameters1.Resid_fric_coeff_layer1) /
                           sl1)) {
        f[Dofw(i, 0)] =
            layer_parameters1.Peak_fric_coeff_layer1 - (sl1 * d[Dofw(i, 0)]);
      } else {
        f[Dofw(i, 0)] = layer_parameters1.Resid_fric_coeff_layer1;
      }

      if (d[Dofw(i, 1)] < ((layer_parameters1.Peak_fric_coeff_layer1 -
                            layer_parameters1.Resid_fric_coeff_layer1) /
                           sl1)) {
        f[Dofw(i, 1)] =
            layer_parameters1.Peak_fric_coeff_layer1 - (sl1 * d[Dofw(i, 1)]);
      } else {
        f[Dofw(i, 1)] = layer_parameters1.Resid_fric_coeff_layer1;
      }

    } else if (id_layers[i] == layer_parameters2.id_layer2) {

      if (d[Dofw(i, 0)] < ((layer_parameters2.Peak_fric_coeff_layer2 -
                            layer_parameters2.Resid_fric_coeff_layer2) /
                           sl2)) {
        f[Dofw(i, 0)] =
            layer_parameters2.Peak_fric_coeff_layer2 - (sl2 * d[Dofw(i, 0)]);
      } else {
        f[Dofw(i, 0)] = layer_parameters2.Resid_fric_coeff_layer2;
      }

      if (d[Dofw(i, 1)] < ((layer_parameters2.Peak_fric_coeff_layer2 -
                            layer_parameters2.Resid_fric_coeff_layer2) /
                           sl1)) {
        f[Dofw(i, 1)] =
            layer_parameters2.Peak_fric_coeff_layer2 - (sl2 * d[Dofw(i, 1)]);
      } else {
        f[Dofw(i, 1)] = layer_parameters2.Resid_fric_coeff_layer2;
      }

    } else if (id_layers[i] == layer_parameters3.id_layer3) {

      if (d[Dofw(i, 0)] < ((layer_parameters3.Peak_fric_coeff_layer3 -
                            layer_parameters3.Resid_fric_coeff_layer3) /
                           sl2)) {
        f[Dofw(i, 0)] =
            layer_parameters3.Peak_fric_coeff_layer3 - (sl3 * d[Dofw(i, 0)]);
      } else {
        f[Dofw(i, 0)] = layer_parameters3.Resid_fric_coeff_layer3;
      }

      if (d[Dofw(i, 1)] < ((layer_parameters3.Peak_fric_coeff_layer3 -
                            layer_parameters3.Resid_fric_coeff_layer3) /
                           sl1)) {
        f[Dofw(i, 1)] =
            layer_parameters3.Peak_fric_coeff_layer3 - (sl2 * d[Dofw(i, 1)]);
      } else {
        f[Dofw(i, 1)] = layer_parameters3.Resid_fric_coeff_layer3;
      }
    }
  }
  return f;
};
}
