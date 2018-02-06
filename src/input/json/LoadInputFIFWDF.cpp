//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 06.02.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include "LoadInputFIFWDF.h"

#include <src/input/json/LoadInputFIFWDF.h>

namespace hfp2d {

using json = nlohmann::json;

//------------------------------------------------------------------------------
// Loading Fluid properties
hfp2d::Fluid loadFluidProperties(json &j_fluid) {
  if (j_fluid.count("Rheology") != 1) {
    std::cout << "No fluid rehology input !";
    il::abort();
  }

  if (j_fluid.count("Density") != 1) {
    std::cout << "No fluid density input !";
    il::abort();
  }
  auto dens = j_fluid["Density"].get<double>();

  if (j_fluid.count("Compressibility") != 1) {
    std::cout << "No fluid Compressibility input !";
    il::abort();
  }
  auto c_f = j_fluid["Compressibility"].get<double>();

  if (j_fluid.count("Rheology") != 1) {
    std::cout << "No fluid rheology input !";
    il::abort();
  }
  // compare Rheology string to Newtonian ?

  if (j_fluid.count("Viscosity") != 1) {
    std::cout << "No fluid viscosity input !";
    il::abort();
  }
  auto visc = j_fluid["Viscosity"].get<double>();

  hfp2d::Fluid the_fluid(dens, c_f, visc);
  return the_fluid;
};
//------------------------------------------------------------------------------
// Loading Solid properties
hfp2d::SolidProperties loadSolidProperties(json &j_rock) {
  if (j_rock.count("Young's modulus") != 1) {
    std::cout << "No Young's modulus input !";
    il::abort();
  }
  auto yme = j_rock["Young's modulus"].get<double>();

  if (j_rock.count("Poisson's ratio") != 1) {
    std::cout << "No Poisson's ratio input !";
    il::abort();
  }
  auto nu = j_rock["Poisson's ratio"].get<double>();

  hfp2d::ElasticProperties my_elas(yme, nu);

  if (j_rock.count("Fracture toughness") != 1) {
    std::cout << "No Fracture toughness input !";
    il::abort();
  }

  if (j_rock.count("Minimum hydraulic width") != 1) {
    std::cout << "No Minimum hydraulic width input !";
    il::abort();
  }

  if (j_rock.count("Leak-off coefficient") != 1) {
    std::cout << "No Leak-off coefficient input !";
    il::abort();
  }

  if (j_rock.count("Increment hydraulic width") != 1) {
    std::cout << "No Increment hydraulic width coefficient input !";
    il::abort();
  }

  if (j_rock.count("Minimum permeability") != 1) {
    std::cout << "No Minimum permeability coefficient input !";
    il::abort();
  }

  if (j_rock.count("Increment permeability") != 1) {
    std::cout << "No Increment permeability coefficient input !";
    il::abort();
  }

  if (j_rock.count("Residual slip") != 1) {
    std::cout << "No Residual slip coefficient input !";
    il::abort();
  }

  if (j_rock.count("Peak friction coefficient") != 1) {
    std::cout << "No Peak friction coefficient coefficient input !";
    il::abort();
  }

  if (j_rock.count("Residual friction coefficient") != 1) {
    std::cout << "No Residual friction coefficient input !";
    il::abort();
  }

  il::int_t nct = j_rock["Fracture toughness"].size();
  il::int_t ncwh = j_rock["Minimum hydraulic width"].size();
  il::int_t ncc = j_rock["Leak-off coefficient"].size();
  il::int_t ncdwh = j_rock["Increment hydraulic width"].size();
  il::int_t nckf = j_rock["Minimum permeability"].size();
  il::int_t ncdkf = j_rock["Increment permeability"].size();
  il::int_t ncd = j_rock["Residual slip"].size();
  il::int_t ncfp = j_rock["Peak friction coefficient"].size();
  il::int_t ncfr = j_rock["Residual friction coefficient"].size();

  IL_EXPECT_FAST((nct == ncwh) && (ncc == nct));

  il::Array<double> tough(nct, 0.), Cl(ncc, 0.), wh_o(ncwh, 0.), kf_o(nckf, 0.),
      f_p(ncfp, 0.), f_r(ncfr, 0.), d_r(ncd, 0.), d_kf(ncdkf, 0.),
      d_wh(ncdwh, 0.);

  for (il::int_t i = 0; i < nct; i++) {
    tough[i] = j_rock["Fracture toughness"][i];
    Cl[i] = j_rock["Leak-off coefficient"][i];
    wh_o[i] = j_rock["Minimum hydraulic width"][i];
    kf_o[i] = j_rock["Minimum permeability"][i];
    f_p[i] = j_rock["Peak friction coefficient"][i];
    f_r[i] = j_rock["Residual friction coefficient"][i];
    d_r[i] = j_rock["Residual slip"][i];
    d_kf[i] = j_rock["Increment permeability"][i];
    d_wh[i] = j_rock["Increment hydraulic width"][i];
  }

  hfp2d::SolidProperties solid(my_elas, tough, wh_o, Cl, kf_o, f_p, f_r, d_r,
                               d_wh, d_kf);

  return solid;
};
//------------------------------------------------------------------------------

hfp2d::InSituConditions loadInSitu(json &j_insitu) {
  if (j_insitu.count("Type") != 1) {
    std::cout << "No type of in-situ conditions in  input !";
    il::abort();
  }

  auto type = j_insitu["Type"].get<std::string>();

  if (type.compare("uniform") != 0) {
    std::cout << "so far only uniform in-situ stress coded up !";
    il::abort();
  }

  auto sxx = j_insitu["Sxx"].get<double>();
  auto syy = j_insitu["Syy"].get<double>();
  auto sxy = j_insitu["Sxy"].get<double>();
  auto p0 = j_insitu["Pore pressure"].get<double>();

  hfp2d::InSituConditions inSitu(sxx, syy, sxy, p0);

  return inSitu;
}
};