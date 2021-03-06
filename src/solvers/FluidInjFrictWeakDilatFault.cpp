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
  // Initialization of variable needed to load the code's arguments
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
  // Instantiate ElasticProperties object
  hfp2d::ElasticProperties ElasticProperties;
  // Instantiate FluidProperties object
  hfp2d::FluidProperties FluidProperties;
  // Instantiate SolidEvolution object
  hfp2d::SolidEvolution SolidEvolution;
  // Instantiate FractureEvolution object
  hfp2d::FractureEvolution FractureEvolution;

  // Load input data from configuration file (*.toml file) and assign them to
  // the previously instantiated objects
  hfp2d::loadInput(input_filename, il::io, MyMesh, ElasticProperties,
                   FluidProperties, SolidEvolution, FractureEvolution,
                   BackgroundLoadingConditions);

  // Pore pressure perturbation at nodal points needed to activate the shear
  // crack (sum of the ambient pore pressure (time t_0) and pore pressure
  // profile at time t_n+1, with n=0.
  // Remember: the pore pressure is always evaluated at time t_n+1
  // (being t_n the current time) as failure must occur!
  // alpha -> initial fault/fracture diffusivity [Lˆ2/T]
  // t_0plus1 -> starting time
  // Dp -> constant overpressure at the middle of the fault
  double alpha = 10;
  double Dp;
  double t_0plus1;
  double time_step;
  double final_time;
  // Import Simulation_Parameters from configuration file (*.toml)
  il::Status status{};
  auto config = il::load<il::MapArray<il::String, il::Dynamic>>(input_filename,
                                                                il::io, status);
  status.abortOnError();
  il::int_t keyFound;
  keyFound = config.search("Simulation_Parameters");

  if (config.found(keyFound) && config.value(keyFound).isMapArray()) {
    const il::MapArray<il::String, il::Dynamic> &SimulationParametersMap =
        config.value(keyFound).asMapArray();

    Dp = findDouble("constant_overpressure", SimulationParametersMap,
                    input_filename);

    t_0plus1 =
        findDouble("initial_time", SimulationParametersMap, input_filename);

    time_step =
        findDouble("time_step", SimulationParametersMap, input_filename);

    final_time =
        findDouble("final_time", SimulationParametersMap, input_filename);

  } else {
    std::cerr << "ERROR: 'Simulation_Parameters' not found in input file "
              << input_filename << std::endl;
    exit(EXIT_FAILURE);
  }

  il::Array<double> press_init_nodes{MyMesh.numberOfNodes(), 0};
  for (il::int_t j = 0; j < press_init_nodes.size(); ++j) {
    press_init_nodes[j] =
        BackgroundLoadingConditions.getAmbientPorePressure(j) +
        (Dp *
         (erfc(fabs((MyMesh.coordinates(j, 0)) / sqrt(alpha * t_0plus1)) / 2)));
  }

  // Matrix to switch from nodal points to collocation points for all DDs (both
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

  // Call SimulationParameters structure and overwrite its default values
  hfp2d::SimulationParameters SimulationParameters;
  SimulationParameters.ehl_max_its = 200;
  SimulationParameters.ehl_tolerance = 1.e-6;
  SimulationParameters.ehl_relaxation = 1.;
  SimulationParameters.frac_front_max_its = 50;
  SimulationParameters.frac_front_tolerance = 1.e-3;

  // Get the elasticity/influence matrix
  il::Array2D<double> kmat{MyMesh.numberDDDofs(), MyMesh.numberDDDofs(), 0};
  kmat = hfp2d::basic_assembly(MyMesh, ElasticProperties,
                               hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

  // Set the source point, i.e the node in the mesh where the fluid is injected
  il::int_t source_point = MyMesh.numberOfElts() / 2;
  hfp2d::Sources Source(source_point);

  /// Initialization of solution object at time t_0, with pressure perturbation
  /// evaluated at t_0plus1
  il::Array<double> init_opening{2 * MyMesh.numberOfElts(), 0.};
  il::Array<double> init_slip{2 * MyMesh.numberOfElts(), 0.};
  double err_frac_position = 2.;
  double err_opening_dd = 2.;
  double err_shear_dd = 2.;
  double err_press = 2.;
  bool expl_impl = false;  // Set true for explicit/implicit integration scheme
  bool damping_term = false;  // Set true for QD simulation
  double shear_modulus = 1.;
  double shear_wave_vel = 1.;
  double damping_coeff = shear_modulus / (2 * shear_wave_vel);
  il::int_t init_iter_front_position = 0;
  il::int_t init_iter_ehls = 0;

  // Pore pressure at time t_0 at collocation points
  auto press_init_coll = il::dot(from_edge_to_coll_press, press_init_nodes);

  // Get the active set of collocation points at time t_0
  // by checking the Mohr-Coulomb criterion
  il::Array<int> init_failed_set_collpoints{0};
  init_failed_set_collpoints.reserve(2 * MyMesh.numberOfElts());
  for (int j = 0, k = 0; j < 2 * MyMesh.numberOfElts(); ++j) {
    if (BackgroundLoadingConditions.getBackgroundShearStress(j) >
        SolidEvolution.getFricCoeff(j) *
            (BackgroundLoadingConditions.getBackgroundNormalStress(j) -
             press_init_coll[j])) {
      init_failed_set_collpoints.resize(k + 1);
      init_failed_set_collpoints[k] = j;
      k = k + 1;
    }
  }

  // Get active set of elements at time t_0
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

  // Call the constructor -> solution object at time t_0
  hfp2d::Solution SolutionAtTn(
      MyMesh, t_0plus1, time_step, init_opening, init_slip, press_init_nodes,
      BackgroundLoadingConditions.getBackgroundNormalStress(),
      BackgroundLoadingConditions.getBackgroundShearStress(),
      init_active_set_elements, init_iter_front_position, init_iter_ehls,
      err_frac_position, err_opening_dd, err_shear_dd, err_press);

  // Initialization of slippage length at time t_0
  double slipp_length_at_Tn;
  if (init_active_set_elements.size() == 0) {
    slipp_length_at_Tn = 0.;
  } else {
    slipp_length_at_Tn = hfp2d::euclidean_distance(
        MyMesh.coordinates(init_active_set_elements[0], 0), 0,
        MyMesh.coordinates(
            init_active_set_elements[init_active_set_elements.size() - 1] + 1,
            0),
        0);
  };

  ////////////////////////        Loop in time          ////////////////////////
  std::string filename;
  double slipp_length_at_Tn_plus1;
  double current_crack_velocity;

  while ((SolutionAtTn.time() < final_time) &&
         (slipp_length_at_Tn <
          euclidean_distance(MyMesh.coordinates(0, 0), 0,
                             MyMesh.coordinates(MyMesh.numberOfNodes() - 1, 0),
                             0)) &&
         (SolutionAtTn.frontIts() < SimulationParameters.frac_front_max_its)) {
    std::cout << " --------------------------------" << std::endl;
    std::cout << " Current time t = " << SolutionAtTn.time() << "\n";

    hfp2d::Solution SolutionAtTnPlus1 = fractFrontPosition(
        kmat, from_edge_to_coll_dds, from_edge_to_coll_dd,
        from_edge_to_coll_press, MyMesh, FluidProperties, SimulationParameters,
        SolidEvolution, FractureEvolution, Source, SolutionAtTn, expl_impl,
        damping_term, damping_coeff);

    // Calculate the new slippage length (i.e at time T_n+1)
    if (SolutionAtTnPlus1.activeElts().size() == 0) {
      slipp_length_at_Tn_plus1 = 0.;
    } else {
      slipp_length_at_Tn_plus1 = hfp2d::euclidean_distance(
          MyMesh.coordinates(SolutionAtTnPlus1.activeElts(0), 0), 0,
          MyMesh.coordinates(SolutionAtTnPlus1.activeElts(
                                 SolutionAtTnPlus1.activeElts().size() - 1) +
                                 1,
                             0),
          0);
    }

    // Calculate the current crack velocity
    if (slipp_length_at_Tn_plus1 - slipp_length_at_Tn == 0) {
      current_crack_velocity = 0.;
    } else {
      current_crack_velocity = (slipp_length_at_Tn_plus1 - slipp_length_at_Tn) /
                               SolutionAtTn.timestep();
    }

    if ((current_crack_velocity > ((2 * 0.02) / SolutionAtTn.timestep())) &&
        (SolutionAtTnPlus1.ehlIts() > SimulationParameters.ehl_max_its / 4)) {
      std::cout
          << "############# Turn into Explicit/Implicit scheme #############"
          << "\n"
          << std::endl;
      expl_impl = true;
    }

    // Write solution at current time T_n to .Json file
    filename = path_output_directory.asCString() +
               std::string{"/Solution_at_time_"} +
               std::to_string(SolutionAtTn.time()) + std::string{".json"};
    SolutionAtTnPlus1.writeToFile(filename);

    // Set new time
    SolutionAtTnPlus1.setTime(SolutionAtTn.time() +
                              SolutionAtTnPlus1.timestep());

    // Update -> old is new
    slipp_length_at_Tn = slipp_length_at_Tn_plus1;
    SolutionAtTn = SolutionAtTnPlus1;
  }
};

///////// FUNCTION FOR FRACTURE FRONT POSITION //////////
Solution fractFrontPosition(
    il::Array2D<double> &elast_matrix, il::Array2D<double> &fetc_dds,
    il::Array2D<double> &fetc_dd, il::Array2D<double> &fetc_press,
    Mesh &theMesh, FluidProperties &FluidProperties,
    SimulationParameters &SimulationParameters, SolidEvolution &SolidEvolution,
    FractureEvolution &FractureEvolution, Sources &Source,
    Solution &SolutionAtTn, bool expl_impl, bool damping_term,
    double damping_coeff) {
  // TODO: insert il_io in the input!

  // Initialization of fracture front loop
  // SolutionAtTn -> solution object at current time Tn
  // SolutionAtTn_k -> current solution object at iteration k of fracture
  // front loop
  // SolutionAtTn_kPlus1 -> solution object at iteration k+1 of fracture front
  // loop
  hfp2d::Solution SolutionAtTn_k = SolutionAtTn;
  hfp2d::Solution SolutionAtTn_kPlus1;
  il::Status status;
  il::Norm norm;
  norm = il::Norm::L1;
  il::int_t cvg_front_posit = 0;  // convergence front position
  il::int_t iter_front_posit = 0;
  il::Array<int> dofs_active_elmnts{};
  dofs_active_elmnts.reserve(theMesh.numberDDDofs());
  il::Array<double> incrm_shearDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openingDD{2 * theMesh.numberOfElts(), 0.};

  while (cvg_front_posit != 1 &&
         iter_front_posit < SimulationParameters.frac_front_max_its) {
    ++iter_front_posit;

    std::cout << "  - Iter for fracture front position = " << iter_front_posit
              << "\n"
              << std::endl;

    // Find the corresponding DOFs of the active elements
    for (int elmnt_i = 0, k = 0; elmnt_i < SolutionAtTn_k.activeElts().size();
         ++elmnt_i) {
      dofs_active_elmnts.resize(k + 4);
      for (int j = 0, l = k; j < theMesh.numberDDDofsPerElt(); ++j, ++l) {
        dofs_active_elmnts[l] =
            theMesh.dofDD(SolutionAtTn_k.activeElts(elmnt_i), j);
      }
      k = k + 4;
    }

    // Call the Reynolds solver for P1 elements
    SolutionAtTn_kPlus1 = hfp2d::reynoldsP1(
        theMesh, elast_matrix, fetc_dds, fetc_dd, fetc_press, SolutionAtTn,
        SolutionAtTn_k, incrm_shearDD, incrm_openingDD, SimulationParameters,
        FluidProperties, SolidEvolution, FractureEvolution, Source,
        dofs_active_elmnts, status, norm, damping_term, damping_coeff);

    // Update active set of elements
    il::Array<int> new_active_set_elements{};
    new_active_set_elements.reserve(theMesh.numberOfElts());
    new_active_set_elements = SolutionAtTn_kPlus1.activeSetElements(
        theMesh, SolutionAtTn_kPlus1, SolidEvolution, fetc_press,
        SolutionAtTn.pressure());

    auto old_active_set_elements = SolutionAtTn_kPlus1.activeElts();
    if (new_active_set_elements.size() == 0) {
      cvg_front_posit = 1;
    } else {
      for (il::int_t i = 0; i < new_active_set_elements.size(); ++i) {
        old_active_set_elements.append(new_active_set_elements[i]);
      }
      std::sort(old_active_set_elements.begin(), old_active_set_elements.end());
      SolutionAtTn_kPlus1.setActiveElts(old_active_set_elements);
    }

    // If the Explcit/Implicit integration scheme is considered, then don't
    // iterate for fracture front position
    if (expl_impl) {
      cvg_front_posit = 1;
    }

    // Set number of iterations of front position to SolutionAtTn_kPlus1
    SolutionAtTn_kPlus1.setFrontPositIters(iter_front_posit);

    // Update -> old is new
    SolutionAtTn_k = SolutionAtTn_kPlus1;
  }

  hfp2d::Solution SolutionAtTnPlus1 = SolutionAtTn_kPlus1;
  return SolutionAtTnPlus1;
};
}
