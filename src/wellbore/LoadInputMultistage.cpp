//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 07.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
#include <string>

#include <il/math.h>

#include <src/wellbore/LoadInputMultistage.h>

namespace hfp2d {

using json = nlohmann::json;

//------------------------------------------------------------------------------
// loading well mesh
hfp2d::WellMesh loadWellMesh(json &j_wmesh) {
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

hfp2d::Sources loadWellSource(json &j_params,hfp2d::WellMesh &the_well){

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
  il::Array<double> hf_vol_rate(n_cl, 0.);  // zero HF rate.

  hfp2d::Sources wellsource(cluster_locs,hf_vol_rate);
  return  wellsource;

}

//------------------------------------------------------------------------------
// loading well parameters.
hfp2d::WellInjection loadWellParameters(json &j_params,
                                        hfp2d::WellMesh &the_well) {

  long cl_p = j_params.count("Clusters MD");  //
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

  IL_EXPECT_FAST((ncp == nct) && (n_cl == ncp) && (ncb == ncp) );

  il::Array<double> perf_coef(ncp, 0.), tort_coef(ncp, 0.), tort_beta(ncp, 0.);

  for (il::int_t i = 0; i < ncp; i++) {
    perf_coef[i] = j_params["Perforations coefficient"][i];
    tort_coef[i] = j_params["Tortuosity coefficient"][i];
    tort_beta[i] = j_params["Tortuosity exponent"][i];
  }

  auto f_model = j_params["Friction model"].get<std::string>();

  il::Array<double> hf_vol_rate(ncp, 0.);  // zero HF rate.

  hfp2d::WellInjection w_inj(pump_rate,  perf_coef,
                             tort_coef, tort_beta);

  return w_inj;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Loading Fluid Properties
hfp2d::Fluid loadFluidProperties(json &j_fluid){

  if (j_fluid.count("Rheology")!=1) {
    std::cout << "No fluid rehology input !";
    il::abort();
  }

  if (j_fluid.count("Density")!=1){
    std::cout << "No fluid density input !";
    il::abort();
  }
  auto dens = j_fluid["Density"].get<double>();

  if (j_fluid.count("Compressibility")!=1){
    std::cout << "No fluid Compressibility input !";
    il::abort();
  }
  auto c_f = j_fluid["Compressibility"].get<double>();

  if (j_fluid.count("Rheology")!=1){
    std::cout << "No fluid rheology input !";
    il::abort();
  }
  // compare Rheology string to Newtonian ?

  if (j_fluid.count("Viscosity")!=1){
    std::cout << "No fluid viscosity input !";
    il::abort();
  }
  auto visc = j_fluid["Viscosity"].get<double>();

  hfp2d::Fluid the_fluid(dens,c_f,visc);
  return the_fluid;
};
//------------------------------------------------------------------------------
// Loading Solid Properties
hfp2d::SolidProperties loadSolidProperties(json &j_rock){

  if (j_rock.count("Young's modulus")!=1) {
    std::cout << "No Young's modulus  input !";
    il::abort();
  }
  auto yme = j_rock["Young's modulus"].get<double>();

  if (j_rock.count("Poisson's ratio")!=1) {
    std::cout << "No Poisson's ratio  input !";
    il::abort();
  }
  auto nu = j_rock["Poisson's ratio"].get<double>();

  hfp2d::ElasticProperties my_elas(yme,nu);

  if (j_rock.count("Fracture toughness")!=1) {
    std::cout << "No Fracture toughness  input !";
    il::abort();
  }

  if (j_rock.count("Minimum hydraulic width")!=1) {
    std::cout << "No Minimum hydraulic width  input !";
    il::abort();
  }

  if (j_rock.count("Leak-off coefficient")!=1) {
    std::cout << "No Leak-off coefficient  input !";
    il::abort();
  }

  il::int_t nct = j_rock["Fracture toughness"].size();
  il::int_t ncwh = j_rock["Minimum hydraulic width"].size();
  il::int_t ncc = j_rock["Leak-off coefficient"].size();

  IL_EXPECT_FAST((nct==ncwh)&&(ncc==nct));


  il::Array<double> tough(nct, 0.), Cl(nct, 0.), wh_o(nct, 0.);

  for (il::int_t i=0;i<nct;i++){
    tough[i]=j_rock["Fracture toughness"][i];
    Cl[i]=j_rock["Leak-off coefficient"][i];
    wh_o[i]=j_rock["Minimum hydraulic width"][i];
  }

  hfp2d::SolidProperties solid(my_elas,tough,wh_o,Cl);

  return solid;

};
//------------------------------------------------------------------------------

hfp2d::InSituConditions loadInSitu(json &j_insitu){

  if (j_insitu.count("Type")!=1) {
    std::cout << "No type of in-situ conditions in  input !";
    il::abort();
  }

  auto type = j_insitu["Type"].get<std::string>();

  if (type.compare("uniform") != 0){
    std::cout << "so far only uniform in-situ stress coded up !";
    il::abort();
  }

  auto sxx = j_insitu["Sxx"].get<double>();
  auto syy = j_insitu["Syy"].get<double>();
  auto sxy = j_insitu["Sxy"].get<double>();


  hfp2d::InSituConditions inSitu(sxx,syy,sxy);

  return inSitu;

}



};
