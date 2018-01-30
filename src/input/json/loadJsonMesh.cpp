//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 28.12.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//


#include <fstream>
#include <src/util/json.hpp>


#include <src/core/Mesh.h>
#include "loadJsonMesh.h"

namespace hfp2d {

using json = nlohmann::json;

hfp2d::Mesh loadJsonMesh(json &j) {
 // std::ifstream input(meshfilename);  // ?const std::string &meshfilename
//  json j;
//  input >> j;

  //  il::Array2D<double> nodes =j["Node coordinates"];
  //
  // check if present.

  int p_present = j.count("Interpolation order");  //
  int n_present = j.count("Node coordinates");     //
  int c_present = j.count("Connectivity");         //
  int m_present = j.count("Material ID");          //
  int f_present = j.count("Fracture ID");          //

  il::int_t inter = 0;  // default
  if (p_present == 1) {
    inter = j["Interpolation order"].get<int>();
  };

  if (n_present != 1) {
    std::cout << " error in Mesh file - no coordinates \n";
    il::abort();
  }

  if (c_present != 1) {
    std::cout << " error in Mesh file - no connectivity \n";
    il::abort();
  }

  il::int_t nn = j["Node coordinates"].size();
  IL_EXPECT_FAST(j["Node coordinates"][0].size() == 2);
  il::Array2D<double> nodes{nn, 2, 0.};

  for (il::int_t i = 0; i < nn; i++) {
    nodes(i, 0) = j["Node coordinates"][i][0];
    nodes(i, 1) = j["Node coordinates"][i][1];
    //    std::cout << "nodes " << i << ": " << nodes(i, 0) << " " << nodes(i,
    //    1)
    //              << "\n";
  }

  il::int_t ne = j["Connectivity"].size();
  IL_EXPECT_FAST(j["Connectivity"][0].size() == 2);

  il::Array2D<il::int_t> conn{ne, 2, 0};
  for (il::int_t i = 0; i < ne; i++) {
    conn(i, 0) = j["Connectivity"][i][0];
    conn(i, 1) = j["Connectivity"][i][1];
    //    std::cout << "connect " << i << ": " << conn(i, 0) << " " << conn(i,
    //    1)
    //              << "\n";
  }

  il::Array<il::int_t> mat{ne,0};

  if (m_present == 1) {  // if a material ID vector is present.
    il::int_t nm = j["Material ID"].size();
    IL_EXPECT_FAST(nm == ne);
    for (il::int_t i = 0; i < nm; i++) {
      mat[i] = j["Material ID"][i];
      //      std::cout << "mat id" << mat[i] <<" \n";
    }

  }

  il::Array<il::int_t> fID{ne,0};

  if (f_present == 1) {  // if a fracture ID vector is present.
    il::int_t nf = j["Fracture ID"].size();
    IL_EXPECT_FAST(nf == ne);
    for (il::int_t i = 0; i < ne; i++) {
      fID[i] = j["Fracture ID"][i];
    }
  }

  hfp2d::Mesh loadedmesh(nodes, conn, mat, fID, inter );

  return loadedmesh;
}
}
