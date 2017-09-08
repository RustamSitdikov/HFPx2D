//
// This file is part of HFPx3D.
//
// Created by D. Nikolski on 1/27/2017.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef INC_3D_BEM_LOAD_MESH_FROM_FILE_H
#define INC_3D_BEM_LOAD_MESH_FROM_FILE_H

#include <cstdio>
#include <il/Array.h>
#include <il/Array2D.h>

namespace hfp2d {

void save_data_to_csv(const il::Array<double> &vector,
                      const il::String &trg_dir, const il::String &of_name) {
  il::String f_path = trg_dir + of_name;
  const char *format = "%.16g\n";
  FILE *of = std::fopen(f_path.c_str(), "w");
  for (int j = 0; j < vector.size(); ++j) {
    std::fprintf(of, format, vector[j]);
  }
  std::fclose(of);
}

void save_data_to_csv(const il::Array2D<double> &matrix,
                      const il::String &trg_dir, const il::String &of_name) {
  il::String f_path = trg_dir + of_name;
  const char *format = "%.16g";
  FILE *of = std::fopen(f_path.c_str(), "w");
  for (int j = 0; j < matrix.size(0); ++j) {
    for (int k = 0; k < matrix.size(1); ++k) {
      double out = matrix(j, k);
      std::fprintf(of, format, out);
      if (k < matrix.size(1) - 1)
        std::fprintf(of, ",");
    }
    std::fprintf(of, "\n");
  }
  std::fclose(of);
}

void export_results(Results_one_timeincrement &SolutionAtTj, double t,
                    const il::String &trg_dir, const il::String &of_name) {

  il::String f_path = trg_dir + of_name;
  const char *format1 = "%.16g %10";
  const char *format6 = "%6.0d";
  const char *format2 = "N. iteration: %i\n\n";
  const char *format3 = "Current time step: %.5g\n\n";
  const char *format4 = "Current time:\n%2.5g\n\n";
  const char *format5 = "Slippage length:\n%2.5g";

  FILE *of = std::fopen(f_path.c_str(), "w");
  std::fprintf(of, format2, SolutionAtTj.iter, 1);
  std::fprintf(of, format3, SolutionAtTj.dt, 1);
  std::fprintf(of, format4, t);
  std::fprintf(of, format5, SolutionAtTj.slippagezone);

  std::fputs("\n\n******* Pressure profile *******\n", of);
  for (int j = 0; j < SolutionAtTj.P.size(); ++j) {
    std::fprintf(of, format1, SolutionAtTj.P[j]);
  }

  std::fprintf(of, "\n\n******* Total slip *******\n");
  for (int j = 0; j < SolutionAtTj.d_tot.size(); ++j) {
    std::fprintf(of, format1, SolutionAtTj.d_tot[j]);
  }

  std::fprintf(of, "\n\n******* Active set of collocation points *******\n");
  for (int j1 = 0; j1 < SolutionAtTj.active_set_collpoints.size(); ++j1) {
    std::fprintf(of, format6, SolutionAtTj.active_set_collpoints[j1]);
  }

  std::fprintf(of, "\n\n******* Friction coefficient *******\n");
  for (int j = 0; j < SolutionAtTj.friction.size(); ++j) {
    std::fprintf(of, format1, SolutionAtTj.friction[j]);
  }

  std::fprintf(of, "\n\n######################################\n\n");

  std::fclose(of);
}
}

#endif // INC_3D_BEM_LOAD_MESH_FROM_FILE_H
