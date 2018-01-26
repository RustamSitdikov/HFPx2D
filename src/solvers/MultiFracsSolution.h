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

  hfp2d::Sources well_outfluxes_ ;

  il::Array<double> dp_entries_ ;  // pressure drop for each cluster.

  il::int_t Its_well_frac_coupling_=0;

  double err_fluxes_=0.;   // shall we also have the residuals here ?

 public:

  //////////////////////////////////////////////////////////////////////////////
  //  CONSTRUCTORS
  //////////////////////////////////////////////////////////////////////////////
  MultiFracsSolution(hfp2d::Solution &frac_s,
                     hfp2d::WellSolution &well_s,
                     hfp2d::Sources &frac_influxes,
                     hfp2d::Sources &well_outfluxes,
                     il::int_t its,
                     il::Array<double> &dp,
                     double err_f) {

    fracs_solution_ = frac_s;
    well_solution_ = well_s;

    IL_EXPECT_FAST(fracs_solution_.time()==well_solution_.time());

    time_=fracs_solution_.time();

    IL_EXPECT_FAST(frac_influxes.SourceElt().size()==well_outfluxes.SourceElt().size());
    frac_influxes_ = frac_influxes;
    well_outfluxes_ =  well_outfluxes;

    il::int_t  nc = frac_influxes.SourceElt().size();
    for (il::int_t i=0; i<nc; i++){
      IL_EXPECT_FAST(frac_influxes.InjectionRate(i)==well_outfluxes.InjectionRate(i));
    }

    IL_EXPECT_FAST(dp.size() == frac_influxes_.InjectionRate().size());


    dp_entries_ = dp;

    err_fluxes_ = err_f;
    Its_well_frac_coupling_ = its;

  };

  ////////////////////////////////////////////////////////////////////////
  //        public interfaces
  ////////////////////////////////////////////////////////////////////////

  hfp2d::Sources fracSources() const { return  frac_influxes_;} ;
  hfp2d::Sources wellSources() const { return  well_outfluxes_;} ;

  hfp2d::Solution fracSolution() const {return fracs_solution_;};

  hfp2d::WellSolution wellSolution() const {return well_solution_;}

  il::Array<il::int_t> wellClusterElt() const { return well_outfluxes_.SourceElt();} ;

  il::Array<il::int_t> fracInletElts() const { return frac_influxes_.SourceElt();} ;

  il::Array<double> fracFluxes() const { return frac_influxes_.InjectionRate();} ;
  double fracFluxes(il::int_t k) const { return frac_influxes_.InjectionRate(k);} ;

  il::Array<double> clusterFluxes() const { return well_outfluxes_.InjectionRate();} ;
  double clusterFluxes(il::int_t k) const { return well_outfluxes_.InjectionRate(k);} ;

  il::Array<double> dpEntries() const { return dp_entries_;} ;
  double dpEntries(il::int_t k) const { return dp_entries_[k];} ;

  double time() const { return time_;};

  //////////////////////////////////////////////////////////////////////////////
  //  METHODS
  //////////////////////////////////////////////////////////////////////////////

  // write to file
  using json = nlohmann::json;

  void writeToFile(std::string &filename) {

    json j_frac = fracs_solution_.createJsonObject() ;

    json j_well = well_solution_.createJsonObject();

    json j_rates = json::array();
    json j_dp = json::array();
    for (il::int_t i=0;i<frac_influxes_.InjectionRate().size();i++){
      j_rates[i]=frac_influxes_.InjectionRate(i);
      j_dp[i]=dp_entries_[i];
    }

    json j_obj= {{"Time", time_},
                 {"Fractures influx",j_rates},
                 {"Clusters Dp entry",j_dp},
                 {"Fractures solution",j_frac},
                 {"Well solution",j_well},
        {"Error well-flow coupling",err_fluxes_},
        {"Iterations well-flow coupling",Its_well_frac_coupling_}
    };

    // write prettified JSON to file
    std::ofstream output(filename);
    output << std::setw(4) << j_obj << std::endl;

  };
};
}
#endif  // HFPX2D_MULTIFRACSSOLUTION_H
