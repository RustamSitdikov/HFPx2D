//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 08.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <algorithm>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>

// Inclusion from the project
#include <il/linear_algebra/dense/norm.h>
#include <src/core/Solution.h>
#include <src/core_dev/FractureEvolution.h>
#include <src/devt/FromEdgeToCol.h>
#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/PlaneStrainInfinite.h>
#include <src/elhsolvers/ReynoldsP1.h>
#include <src/input/LoadArguments.h>
#include <src/input/LoadInput.h>

namespace hfp2d {

void fluidInjFrictWeakDilatFault(int argc, char const *argv[]) {
  // Initialization
  bool check_input = false;
  bool check_output = false;
  bool check_restart = false;
  il::String input_filename;
  il::String restart_filename;
  il::String path_output_directory;

  // Load arguments & create the output directory
  hfp2d::loadArguments(argc, argv, il::io, check_input, check_output,
                       check_restart, input_filename, restart_filename,
                       path_output_directory);

  // Intantiate mesh object
  hfp2d::Mesh MyMesh;
  // Instantiate InSituStress object
  hfp2d::InSituStress BackgroundLoadingConditions;
  // Instantiate elasticProperties object
  hfp2d::ElasticProperties ElasticProperties;
  // Instantiate FluidProperties object
  hfp2d::FluidProperties FluidProperties;
  // Instantiate SolidEvolution object
  hfp2d::SolidEvolution SolidEvolution;
  // Instantiate FractureEvolution object
  hfp2d::FractureEvolution FractureEvolution;

  // Load input data from .toml file
  hfp2d::loadInput(input_filename, il::io, MyMesh, ElasticProperties,
                   FluidProperties, SolidEvolution, FractureEvolution,
                   BackgroundLoadingConditions);

  // Pore pressure profile at nodal points needed to activate the shear crack
  // (sum of the ambient pore pressure and pore pressure profile at time t_n+1,
  // where n=0
  // Remember that pore pressure is always evaluated at time t_n+1 as failure
  // must occur in order to activate a shear crack)
  // alpha -> initial fault/fracture diffusivity [Lˆ2/T]
  // t_0plus -> initial time (starting time)
  double alpha = 10;
  double Dp = 0.5;
  double t_0plus = 0.00003;
  double init_dt = 0.00003;
  double time_max = 0.001;
  il::Array<double> press_init_nodes{MyMesh.numberOfNodes(), 0};
  for (il::int_t j = 0; j < press_init_nodes.size(); ++j) {
    press_init_nodes[j] =
        BackgroundLoadingConditions.getAmbientPorePressure(j) +
        (Dp *
         (erfc(fabs((MyMesh.coordinates(j, 0)) / sqrt(alpha * t_0plus)) / 2)));
  }

  // Matrix to switch from nodal points to collocation points for DDs (both
  // opening and shear)
  il::Array2D<double> from_edge_to_coll_dds{MyMesh.numberDDDofs(),
                                            MyMesh.numberDDDofs(), 0};
  from_edge_to_coll_dds = hfp2d::from_edge_to_col_dg_full2d(MyMesh);

  // Matrix to switch from nodal points to collocation points for one DD (either
  // opening or shear)
  il::Array2D<double> from_edge_to_coll_dd{MyMesh.numberDDDofs(),
                                           MyMesh.numberDDDofs(), 0};
  from_edge_to_coll_dd = hfp2d::from_edge_to_col_dg(MyMesh);

  // Matrix to switch from nodal points to collocation points for pressure
  il::Array2D<double> from_edge_to_coll_press{MyMesh.numberDDDofs(),
                                              MyMesh.numberOfNodes(), 0};
  from_edge_to_coll_press = hfp2d::from_edge_to_col_cg(MyMesh);

  // Call SimulationParameters structure with its default values
  hfp2d::SimulationParameters SimulationParameters;

  // Get the elasticity/influence matrix
  il::Array2D<double> kmat{MyMesh.numberDDDofs(), MyMesh.numberDDDofs(), 0};
  kmat = hfp2d::basic_assembly(MyMesh, ElasticProperties,
                               hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

  // Set the source point, i.e the node in the mesh where the fluid is injected
  il::int_t source_point = MyMesh.numberOfElts() / 2;
  hfp2d::Sources Source(source_point);

  /// Call constructor of class Solution -> Initialization of SolutionAtTn
  /// object
  il::Array<double> init_opening{2 * MyMesh.numberOfElts(), 0.};
  il::Array<double> init_slip{2 * MyMesh.numberOfElts(), 0.};
  double err_frac_position = 2.;
  double err_opening_dd = 2.;
  double err_shear_dd = 2.;
  double err_press = 2.;
  il::int_t init_iter_front_position = 0;
  il::int_t init_iter_ehls = 0;

  il::Array<double> press_init_coll{2 * MyMesh.numberOfElts(), 0};
  auto p_init_coll = il::dot(from_edge_to_coll_press, press_init_nodes);
  for (il::int_t i = 0, k = 1; i < press_init_coll.size(); ++i, k = k + 2) {
    press_init_coll[i] = p_init_coll[k];
  }

  // Get the active set of collocation points by checking the MC criterion
  il::Array<int> init_failed_set_collpoints{0};
  init_failed_set_collpoints.reserve(2 * MyMesh.numberOfElts());
  for (int j = 0, k = 0; j < 2 * MyMesh.numberOfElts(); ++j) {
    if (BackgroundLoadingConditions.getBackgroundShearStress(j) >=
        SolidEvolution.getFricCoeff(j) *
            (BackgroundLoadingConditions.getBackgroundNormalStress(j) -
             press_init_coll[j])) {
      init_failed_set_collpoints.resize(k + 1);
      init_failed_set_collpoints[k] = j;
      k = k + 1;
    }
  }

  // Get initial active set of elements
  il::Array2D<int> dof_single_dd{MyMesh.numberOfElts(),
                                 (MyMesh.interpolationOrder() + 1), 0};
  for (int i = 0; i < MyMesh.numberOfElts(); i++) {
    for (int j = 0; j < 1 * (MyMesh.interpolationOrder() + 1); j++) {
      dof_single_dd(i, j) = i * 1 * (MyMesh.interpolationOrder() + 1) + j;
    }
  }

  il::Array<int> init_set_elements{0};
  init_set_elements.reserve(2 * MyMesh.numberOfElts());
  for (int l = 0, k = 0; l < init_failed_set_collpoints.size(); ++l, ++k) {
    init_set_elements.resize(k + 1);
    init_set_elements[l] =
        hfp2d::find_2d_integer(dof_single_dd, init_failed_set_collpoints[l])[0];
  }
  auto init_set_elmnts = hfp2d::delete_duplicates_integer(init_set_elements);

  il::Array<int> init_active_set_elements{init_set_elmnts.size()};

  if (init_failed_set_collpoints.size() ==
      2 * init_active_set_elements.size()) {
    init_active_set_elements = init_set_elmnts;
  } else {
    init_active_set_elements.resize(init_set_elmnts.size() - 2);
    for (int i = 0, k = 1; i < init_active_set_elements.size(); ++i, ++k) {
      init_active_set_elements[i] = init_set_elmnts[k];
    }
  }

  // Initialization of slippage length
  double slipp_length;
  if (init_active_set_elements.size() == 0) {
    slipp_length = 0.;

  } else {
    slipp_length = hfp2d::euclidean_distance(
        MyMesh.coordinates(init_active_set_elements[0], 0), 0,
        MyMesh.coordinates(
            init_active_set_elements[init_active_set_elements.size() - 1] + 1,
            0),
        0);
  };

  // Constructor
  hfp2d::Solution SolutionAtTn(
      MyMesh, t_0plus, init_dt, init_opening, init_slip, press_init_nodes,
      BackgroundLoadingConditions.getBackgroundNormalStress(),
      BackgroundLoadingConditions.getBackgroundShearStress(),
      init_active_set_elements, init_iter_front_position, init_iter_ehls,
      err_frac_position, err_opening_dd, err_shear_dd, err_press);

  /// Loop in time
  std::string filename;
  while (SolutionAtTn.time() <= time_max &&
         slipp_length <= euclidean_distance(
                             MyMesh.coordinates(0, 0), 0,
                             MyMesh.coordinates(MyMesh.numberOfNodes() - 1, 0),
                             0) &&
         SolutionAtTn.frontIts() <= SimulationParameters.frac_front_max_its) {
    std::cout << "******** Current time ******* "
              << "t = " << SolutionAtTn.time() << "\n";

    SolutionAtTn = fractFrontPosition(
        kmat, from_edge_to_coll_dds, from_edge_to_coll_dd,
        from_edge_to_coll_press, MyMesh, FluidProperties, SimulationParameters,
        SolidEvolution, FractureEvolution, Source, SolutionAtTn);

    if (SolutionAtTn.activeElts().size() == 0) {
      slipp_length = 0.;
    } else {
      slipp_length = hfp2d::euclidean_distance(
          MyMesh.coordinates(SolutionAtTn.activeElts(0), 0), 0,
          MyMesh.coordinates(
              SolutionAtTn.activeElts(init_active_set_elements.size() - 1) + 1,
              0),
          0);
    }

    filename = std::string{"/Users/federicociardo/Desktop/Solution_at_time_"} +
               std::to_string(SolutionAtTn.time()) + std::string{".json"};
    SolutionAtTn.writeToFile(filename);
  }
};

///////// FUNCTION FOR FRACTURE FRONT POSITION //////////
Solution fractFrontPosition(il::Array2D<double> &elast_matrix,
                            il::Array2D<double> &fetc_dds,
                            il::Array2D<double> &fetc_dd,
                            il::Array2D<double> &fetc_press, Mesh &theMesh,
                            FluidProperties &FluidProperties,
                            SimulationParameters &SimulationParameters,
                            SolidEvolution &SolidEvolution,
                            FractureEvolution &FractureEvolution,
                            Sources &Source, Solution &SolutionAtTn) {
  // TODO: insert il_io in the input!

  // Initialization of fracture front loop
  il::Status status;
  il::Norm norm;
  norm = il::Norm::L1;
  //  Solution SolutionAtnPlusOne;
  il::int_t cvg_front_posit = 0;
  il::int_t iter_front_posit = 1;
  il::Array<int> dof_active_elmnts{};
  dof_active_elmnts.reserve(theMesh.numberDDDofs());

  il::Array<double> tau_old{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t n = 0; n < tau_old.size(); ++n) {
    tau_old[n] = SolutionAtTn.tau(n);
  }

  il::Array<double> sigmaN_old{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t n = 0; n < tau_old.size(); ++n) {
    sigmaN_old[n] = SolutionAtTn.sigmaN(n);
  }

  il::Array<double> shearDD_old{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t n = 0; n < shearDD_old.size(); ++n) {
    shearDD_old[n] = SolutionAtTn.shearDD(n);
  }

  il::Array<double> openingDD_old{2 * theMesh.numberOfElts(), 0.};
  for (il::int_t n = 0; n < openingDD_old.size(); ++n) {
    openingDD_old[n] = SolutionAtTn.openingDD(n);
  }

  il::Array<double> press_old{theMesh.numberOfNodes(), 0.};
  for (il::int_t i = 0; i < press_old.size(); ++i) {
    press_old[i] = SolutionAtTn.pressure(i);
  }

  il::Array<double> incrm_shearDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openingDD{2 * theMesh.numberOfElts(), 0.};

  while (cvg_front_posit != 1 &&
         iter_front_posit <= SimulationParameters.frac_front_max_its) {
    std::cout << "Iter for fracture front position = " << iter_front_posit
              << "\n"
              << std::endl;

    // Find the corresponding DOFs of the active elements
    for (int elmnt_i = 0, k = 0; elmnt_i < SolutionAtTn.activeElts().size();
         ++elmnt_i) {
      dof_active_elmnts.resize(k + 4);
      for (int j = 0, l = k; j < theMesh.numberDDDofsPerElt(); ++j, ++l) {
        dof_active_elmnts[l] =
            theMesh.dofDD(SolutionAtTn.activeElts(elmnt_i), j);
      }
      k = k + 4;
    }

    SolutionAtTn = hfp2d::reynoldsP1(
        theMesh, elast_matrix, fetc_dds, fetc_dd, fetc_press, SolutionAtTn,
        tau_old, sigmaN_old, shearDD_old, openingDD_old, press_old,
        incrm_shearDD, incrm_openingDD, SimulationParameters, FluidProperties,
        SolidEvolution, FractureEvolution, Source, dof_active_elmnts, status,
        norm);

    // Update active set of elements
    il::Array<int> new_active_set_elements{};
    new_active_set_elements.reserve(theMesh.numberOfElts());
    new_active_set_elements = SolutionAtTn.activeSetElements(
        theMesh, SolutionAtTn, SolidEvolution, fetc_press, press_old);

    auto old_active_set_elements = SolutionAtTn.activeElts();

    if (new_active_set_elements.size() == 0) {
      cvg_front_posit = 1;

    } else {
      for (il::int_t i = 0; i < new_active_set_elements.size(); ++i) {
        old_active_set_elements.append(new_active_set_elements[i]);
      }

      std::sort(old_active_set_elements.begin(), old_active_set_elements.end());
      SolutionAtTn.setActiveElts(old_active_set_elements);
    }

    ++iter_front_posit;
  }

  SolutionAtTn.setFrontPositIters(iter_front_posit - 1);

  return SolutionAtTn;
};
}
