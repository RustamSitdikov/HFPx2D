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
#include <src/input/json/LoadInputMultistage.h>

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
    std::cout << "No parameters in input file ";
    il::abort();
  }
  json j_params = j["Model parameters"];

  hfp2d::WellMesh the_well= loadWellMesh(j_wmesh);
  hfp2d::WellInjection w_inj= loadWellParameters(j_params, the_well);
  hfp2d::Sources w_outflow = loadWellSource(j_params,the_well);

  if (j_params.count("Fluid properties")!=1){
    std::cout << "No fluid properties input in  model parameters";
    il::abort();
  }

  json j_fluid = j_params["Fluid properties"];

  hfp2d::Fluid water=loadFluidProperties(j_fluid);

  double dt = 0.1;

  // create initial solution
  hfp2d::WellSolution iniSol(the_well, w_inj, water);
  hfp2d::SimulationParameters WellFlowParam;

  WellFlowParam.ehl_relaxation = 1.0;
  WellFlowParam.ehl_tolerance = 1.e-6;

  // one step.

  hfp2d::WellSolution SolN = iniSol;


  WellSolution SolN_1 = wellFlowSolverP0(SolN,
                                         the_well,
                                         w_inj,
                                         w_outflow,
                                         ffChurchill,
                                         dt,
                                         WellFlowParam,
                                         false,
                                         water);

  std::string resfilename = "../Debug/WellTest_results.json";

  SolN_1.writeToFile(resfilename);

  return 0;
}
}
