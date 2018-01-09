//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 05.01.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#ifndef HFPX2D_MULTIFRACSSOLUTION_H
#define HFPX2D_MULTIFRACSSOLUTION_H

#include <fstream>

#include <src/core/Solution.h>
#include <src/core/Sources.h>
#include <src/wellbore/WellSolution.h>
#include <src/util/json.hpp>

namespace hfp2d {

class MultiFracsSolution {
 private:
  hfp2d::Solution fracs_solution_;

  hfp2d::WellSolution well_solution_;

  double time_;

  hfp2d::Sources frac_influxes_ ;

  // il::Array<double> dp_entries_ ;  // pressure drop for each cluster.

  il::int_t Its_well_frac_coupling_=0;
  double err_fluxes_=0.;

 public:

  //////////////////////////////////////////////////////////////////////////////
  //  CONSTRUCTORS
  //////////////////////////////////////////////////////////////////////////////
  MultiFracsSolution(Solution &frac_s,WellSolution &well_s,
                     Sources &rates_s,double err_f,il::int_t its) {

    fracs_solution_ = frac_s;
    well_solution_ = well_s;

    IL_EXPECT_FAST(fracs_solution_.time()==well_solution_.time());
    time_=fracs_solution_.time();

    frac_influxes_ = rates_s;

    err_fluxes_=err_f;
    Its_well_frac_coupling_=its;

  };

  //////////////////////////////////////////////////////////////////////////////
  //  METHODS
  //////////////////////////////////////////////////////////////////////////////

  // write to file
  using json = nlohmann::json;

  void writeToFile(std::string &filename) {

    json j_frac = fracs_solution_.createJsonObject() ;

    json j_well = well_solution_.createJsonObject();

    json j_rates = json::array();
    for (il::int_t i=0;i<frac_influxes_.InjectionRate().size();i++){
      j_rates=frac_influxes_.InjectionRate(i);
    }

    json j_obj= {{"Time", time_},
                 {"Fractures influx",j_rates},
                 {"Fractures solution",j_frac},
                 {"Well solution",j_well},
        {"Error well-flow coupling",err_fluxes_},
        {"Its well-flow coupling",Its_well_frac_coupling_}
    };


    // write prettified JSON to file
    std::ofstream output(filename);
    output << std::setw(4) << j_obj << std::endl;

  };
};
}
#endif  // HFPX2D_MULTIFRACSSOLUTION_H
