//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 04.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include <fstream>

#include <src/wellbore/SimpleWellFlowSolverBenchmark.h>
#include <src/wellbore/WellFlowP0.h>
#include <src/wellbore/LoadInputMultistage.h>
#include <src/util/json.hpp>

namespace hfp2d {

using json = nlohmann::json;

int WellboreFlowBenchmark() {


  std::string wellfilename = "../Debug/WellTest.json";

  std::ifstream input(wellfilename);  // ?
  json j;
  input >> j;

  if ((j.count("Wellbore mesh") != 1)) {
    std::cout << "No wellbore mesh in json input file ";
    il::abort();
  }
  json j_wmesh = j["Wellbore mesh"];

  if ((j.count("Model parameters") != 1)) {
    std::cout << "No wellbore parameters in input file ";
    il::abort();
  }
  json j_params = j["Model parameters"];

  hfp2d::WellMesh the_well=LoadWellMesh(j_wmesh);
  hfp2d::WellInjection w_inj=LoadWellParameters(j_params,the_well);

  double dt = 0.1;

  // create initial solution
  hfp2d::Fluid water(1.e3, 1.e-3, 5.e-10);

  hfp2d::WellSolution iniSol(the_well, w_inj, water, dt);

  hfp2d::SimulationParameters WellFlowParam;

  WellFlowParam.ehl_relaxation = 1.0;
  WellFlowParam.ehl_tolerance = 1.e-6;
  // one step.

  hfp2d::WellSolution SolN = iniSol;

  hfp2d::WellSolution SolN_1 = wellFlowSolverP0(
      SolN, the_well, w_inj, water, ffChurchill, dt, WellFlowParam, false);

  std::string resfilename = "../Debug/WellTest_results.json";

  SolN_1.writeToFile(resfilename);


  return 0;
}
}
