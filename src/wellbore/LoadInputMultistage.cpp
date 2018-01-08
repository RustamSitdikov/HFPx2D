//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 07.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include <il/math.h>

#include "LoadInputMultistage.h"

#include <src/wellbore/WellInjection.h>
#include <src/wellbore/WellMesh.h>
#include <src/util/json.hpp>

namespace hfp2d {

using json = nlohmann::json;

//------------------------------------------------------------------------------
hfp2d::WellMesh LoadWellMesh(json &j_wmesh) {
  // load the well mesh and create the welll mesh object from the json
  long md_present = j_wmesh.count("MD");               //
  long tvd_present = j_wmesh.count("TVD");             //
  long pipe_present = j_wmesh.count("Pipe diameter");  //
  long r_present = j_wmesh.count("Pipe roughness");    //
  long wa_present = j_wmesh.count("Well azimuth");

  if (md_present != 1) {
    std::cout << " error in wellMesh file - no MD \n";
    il::abort();
  }

  if (tvd_present != 1) {
    std::cout << " error in wellMesh file - no TVD \n";
    il::abort();
  }

  if (pipe_present != 1) {
    std::cout << " error in wellMesh file - no pipe diameters \n";
    il::abort();
  }

  double pipe_r = 0.;  // default  -> currently we only enter a
  if (r_present == 1) {
    pipe_r = j_wmesh["Pipe roughness"].get<double>();
  };

  il::int_t nn = j_wmesh["MD"].size();
  il::int_t nt = j_wmesh["TVD"].size();
  IL_EXPECT_FAST(nn == nt);

  il::Array<double> md(nn, 0.), tvd(nt, 0);
  for (il::int_t i = 0; i < nn; i++) {
    md[i] = j_wmesh["MD"][i];
    tvd[i] = j_wmesh["TVD"][i];
  }

  il::int_t ne = j_wmesh["Pipe diameter"].size();
  IL_EXPECT_FAST(ne == nt - 1);

  il::Array<double> pipediam(ne, 0.), piperough(ne, 0.);
  for (il::int_t i = 0; i < ne; i++) {
    pipediam[i] = j_wmesh["Pipe diameter"][i];
    piperough[i] = pipe_r;
  }

  double az = 0.;
  if (wa_present == 1) {
    az = j_wmesh["Well azimuth"].get<double>();
  }

  hfp2d::WellMesh the_well(md, tvd, pipediam, az, piperough);
  return the_well;
}
//------------------------------------------------------------------------------

hfp2d::WellInjection LoadWellParameters(json &j_params,
                                        hfp2d::WellMesh &the_well) {

  long cl_p = j_params.count("Clusters MD");  //
  // int tvd_present = j_params.count("Plug MD");         //
  long p_present = j_params.count("Perforations coefficient");  //
  long t_present = j_params.count("Tortuosity coefficient");    //
  long t_b_present = j_params.count("Tortuosity exponent");

  if (cl_p != 1) {
    std::cout << " error in input file - no cluster MD \n";
    il::abort();
  }

  il::int_t n_cl = j_params["Clusters MD"].size();
  il::Array<double> cl_mds(n_cl, 0.);

  for (il::int_t i = 0; i < n_cl; i++) {
    cl_mds[i] = j_params["Clusters MD"][i];
  }

  il::Array<il::int_t> cluster_locs(n_cl, 0);
  // localize the element where the cluster is located
  il::int_t n_elts = the_well.md().size() - 1;
  for (il::int_t i = 0; i < n_cl; i++) {
    for (il::int_t e = 0; e < n_elts; e++) {
      if ((cl_mds[i] < the_well.md(the_well.connectivity(e, 1))) &&
          (cl_mds[i] >= the_well.md(the_well.connectivity(e, 0)))) {
        cluster_locs[i] = e;
        break;
      }
    }
  }

  auto pump_rate = j_params["Pump rate"].get<double>();

  // for now we do not have the plug md, we assume that the end of the well mesh
  // correspond to the plug location
  //  auto plug_md = j_params["Plug MD"].get<double>();
  //  il::int_t plug_loc = ne;
  //  for (il::int_t e = 0; e < the_well.numberOfElts(); e++) {
  //    if ((plug_md < md[the_well.connectivity(e, 1)]) &&
  //        (plug_md >= md[the_well.connectivity(e, 0)])) {
  //      plug_loc = e;
  //      break;
  //    }
  //  }

  if (p_present != 1) {
    std::cout << "error in input file - no perf coef \n";
    il::abort();
  }
  if (t_present != 1) {
    std::cout << " error in input file - no tortuosity coef \n";
    il::abort();
  }

  if (t_b_present != 1) {
    std::cout << " error in input file - no tortuosity exponent \n";
    il::abort();
  }

  il::int_t ncp = j_params["Perforations coefficient"].size();
  il::int_t nct = j_params["Tortuosity coefficient"].size();
  il::int_t ncb = j_params["Tortuosity exponent"].size();

  IL_EXPECT_FAST((ncp == nct) && (n_cl == ncp) && (ncb == ncp));

  il::Array<double> perf_coef(ncp, 0.), tort_coef(ncp, 0.), tort_beta(ncp, 0.);

  for (il::int_t i = 0; i < ncp; i++) {
    perf_coef[i] = j_params["Perforations coefficient"][i];
    tort_coef[i] = j_params["Tortuosity coefficient"][i];
    tort_beta[i] = j_params["Tortuosity exponent"][i];
  }

  auto f_model = j_params["Friction model"].get<std::string>();

  il::Array<double> hf_vol_rate(ncp, 0.);  // zero HF rate.

  hfp2d::WellInjection w_inj(pump_rate, cluster_locs, hf_vol_rate, perf_coef,
                             tort_coef, tort_beta);

  return w_inj;
}
};
