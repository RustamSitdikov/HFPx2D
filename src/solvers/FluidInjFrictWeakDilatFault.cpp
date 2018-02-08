//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 08.11.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from standard library
#include <unistd.h>
#include <algorithm>
#include <fstream>
#include <iostream>

// Inclusion from Inside Loop library
#include <il/Array.h>
#include <il/Array2D.h>
#include <il/String.h>

// Inclusion from the project
#include <il/linear_algebra/dense/norm.h>
#include <src/core/Fluid.h>
#include <src/core/SimulationParameters.h>
#include <src/core/Solution.h>
#include <src/core_dev/FractureEvolution.h>
#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/PlaneStrainInfinite.h>
#include <src/input/json/LoadInputFIFWDF.h>
#include <src/input/json/loadJsonMesh.h>
#include <src/solvers/FluidInjFrictWeakDilatFault.h>
#include <src/util/FromEdgeToCol.h>

namespace hfp2d {

////////////////////////////////////////////////////////////////////////////////
using json = nlohmann::json;

void fluidInjFrictWeakDilatFault(int argc, char const *argv[]) {
  /////// LOAD PROGRAM ARGUMENT ///////
  // Initial check
  if (argc != 1) {
    std::cerr << "ERROR -> Usage of the program is: input_file" << std::endl;
    std::cerr << "-- Press ENTER to exit...";
    std::cin.get();
    exit(EXIT_FAILURE);
  }

  // Instantiate an object from stringstream class
  std::stringstream arg_stream;

  // Insert argument in the object
  arg_stream << argv[0];

  // Extract argument from arg_stream object
  std::string program_argument;
  arg_stream >> program_argument;

  il::int_t status_of_access;
  auto config_filename = il::String(
      il::StringType::Ascii, program_argument.c_str(), program_argument.size());

  status_of_access = access(config_filename.asCString(), F_OK | R_OK);

  if (status_of_access != 0) {
    std::cerr << "Error: impossible to access the input file "
              << config_filename << std::endl;
    std::cerr << strerror(errno) << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ifstream input(program_argument);
  json js;
  input >> js;

  // Fracture mesh
  if ((js.count("Fractures mesh") != 1)) {
    std::cout << "No fracture mesh in json input file ";
    il::abort();
  }
  json j_fmesh = js["Fractures mesh"];
  hfp2d::Mesh frac_mesh = loadJsonMesh(j_fmesh);

  // Model parameters
  //  - Fluid properties
  //  - Rock properties
  //  - In-situ conditions
  //  - Constant overpressure
  //  - Injection rate
  if ((js.count("Model parameters") != 1)) {
    std::cout << "No parameters in input file ";
    il::abort();
  }
  json j_params = js["Model parameters"];
  hfp2d::Fluid fracfluid = loadFluidProperties(j_params["Fluid properties"]);

  if (j_params.count("Rock properties") != 1) {
    std::cout << "No rock properties input in model parameters";
    il::abort();
  }
  json j_rock = j_params["Rock properties"];
  hfp2d::SolidProperties rock = loadSolidProperties(j_rock);

  if (j_rock.count("Dilatancy effect on sigmaN") != 1) {
    std::cout << "No Dilatancy effect on sigmaN in rock parameters";
    il::abort();
  }
  il::int_t ndilplast = j_rock["Dilatancy effect on sigmaN"].size();
  il::Array<double> dilat_plast{ndilplast, 0.};
  for (il::int_t i = 0; i < ndilplast; i++) {
    dilat_plast[i] = j_rock["Dilatancy effect on sigmaN"][i];
  }

  if (j_params.count("In-situ conditions") != 1) {
    std::cout << "No In-situ conditions input in  model parameters";
    il::abort();
  }
  json j_insitu = j_params["In-situ conditions"];
  hfp2d::InSituConditions inSituStress = loadInSitu(j_insitu);

  if (j_params.count("Constant overpressure") != 1) {
    std::cout << "No Constant overpressure in input !";
    il::abort();
  }
  double const_overpress = j_params["Constant overpressure"];

  if (j_params.count("Injection rate") != 1) {
    std::cout << "No Injection rate in input !";
    il::abort();
  }
  double inj_rate = j_params["Injection rate"];

  // Simulation parameters
  if ((js.count("Simulation parameters") != 1)) {
    std::cout << "No Simulation parameters in input file ";
    il::abort();
  }
  json j_simul = js["Simulation parameters"];

  if (j_simul.count("Maximum time") != 1) {
    std::cout << "No Maximum time in input !";
    il::abort();
  }
  double t_max = j_simul["Maximum time"];

  if (j_simul.count("Initial time") != 1) {
    std::cout << "No Initial time in input !";
    il::abort();
  }
  double t_ini = j_simul["Initial time"];

  if (j_simul.count("Adaptive time stepping") != 1) {
    std::cout << "No Adaptive time stepping in input !";
    il::abort();
  }
  bool adapt_dt = j_simul["Adaptive time stepping"];

  if (j_simul.count("Initial time step") != 1) {
    std::cout << "No Initial time step in input !";
    il::abort();
  }
  double init_dt = j_simul["Initial time step"];

  if (j_simul.count("Maximum time step") != 1) {
    std::cout << "No Maximum time step in input !";
    il::abort();
  }
  double dt_max = j_simul["Maximum time step"];

  if (j_simul.count("Explicit/Implicit scheme") != 1) {
    std::cout << "No Explicit/Implicit scheme in input !";
    il::abort();
  }
  bool expl_impl = j_simul["Explicit/Implicit scheme"];

  if (j_simul.count("Quasi-Dynamic formulation") != 1) {
    std::cout << "No Quasi-Dynamic formulation in input !";
    il::abort();
  }
  bool QD = j_simul["Quasi-Dynamic formulation"];

  if (j_simul.count("Fracture front tolerance") != 1) {
    std::cout << "No Fracture front tolerance in input !";
    il::abort();
  }
  double frac_front_tolerance = j_simul["Fracture front tolerance"];

  // Output directory
  std::string output_dir;
  if ((js.count("Results files name core") == 1)) {
    output_dir = js["Results files name core"].get<std::string>();
  }

  // Instantiation of SolidEvolution & FractureEvolution object
  il::Array<double> peak_fric_coeff{2 * frac_mesh.numberOfElts(),
                                    rock.fpeak(0)};
  il::Array<double> res_fric_coeff{2 * frac_mesh.numberOfElts(), rock.fres(0)};
  il::Array<double> res_slip{2 * frac_mesh.numberOfElts(), rock.resid_slip(0)};
  hfp2d::SolidEvolution SolidEvolution(peak_fric_coeff, res_fric_coeff,
                                       res_slip);

  il::Array<double> init_permeab{2 * frac_mesh.numberOfElts(), rock.Kf_O(0)};
  il::Array<double> incr_permeab{2 * frac_mesh.numberOfElts(),
                                 rock.incrm_kf(0)};
  il::Array<double> init_hydr_width{2 * frac_mesh.numberOfElts(), rock.Wh_O(0)};
  il::Array<double> incr_hydr_width{2 * frac_mesh.numberOfElts(),
                                    rock.incrm_wh(0)};
  hfp2d::FractureEvolution FractureEvolution(
      init_permeab, incr_permeab, res_slip, init_hydr_width, incr_hydr_width);

  /////// INITIALIZATION OF OTHER VARIABLES ///////
  // Matrix to switch from nodal points to collocation points for all DDs
  // (both opening and shear)
  il::Array2D<double> from_edge_to_coll_dds{frac_mesh.numberDDDofs(),
                                            frac_mesh.numberDDDofs(), 0};
  from_edge_to_coll_dds = hfp2d::from_edge_to_col_dg_full2d(frac_mesh);

  // Matrix to switch from nodal points to collocation points for one DD
  // (either opening or shear)
  il::Array2D<double> from_edge_to_coll_dd{frac_mesh.numberDDDofs(),
                                           frac_mesh.numberDDDofs(), 0};
  from_edge_to_coll_dd = hfp2d::from_edge_to_col_dg(frac_mesh);

  // Matrix to switch from nodal points to collocation points for pressure
  il::Array2D<double> from_edge_to_coll_press{frac_mesh.numberDDDofs(),
                                              frac_mesh.numberOfNodes(), 0};
  from_edge_to_coll_press = hfp2d::from_edge_to_col_cg(frac_mesh);

  // Dof handle for single DD (either shear DD or opening DD)
  il::Array2D<il::int_t> dof_single_dd{frac_mesh.numberOfElts(),
                                       (frac_mesh.interpolationOrder() + 1), 0};
  for (il::int_t i = 0; i < frac_mesh.numberOfElts(); i++) {
    for (il::int_t j = 0; j < 1 * (frac_mesh.interpolationOrder() + 1); j++) {
      dof_single_dd(i, j) = i * 1 * (frac_mesh.interpolationOrder() + 1) + j;
    }
  }

  // Call SimulationParameters structure and overwrite its default values
  hfp2d::SimulationParameters SimulationParameters;
  SimulationParameters.ehl_max_its = 200;
  SimulationParameters.ehl_tolerance = 1.e-6;
  SimulationParameters.ehl_relaxation = 0.7;
  SimulationParameters.frac_front_max_its = 50;
  SimulationParameters.frac_front_tolerance = frac_front_tolerance;

  // Get the elasticity/influence matrix
  hfp2d::ElasticProperties elast_properties = rock.ElasticProperties();

  il::Array2D<double> kmat{frac_mesh.numberDDDofs(), frac_mesh.numberDDDofs(),
                           0};
  kmat = hfp2d::basic_assembly(frac_mesh, elast_properties,
                               hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

  // Set the source point, i.e the node in the mesh where the fluid is injected
  // It's an array as multiple injection points may exist
  il::Array<il::int_t> source_point{1, frac_mesh.numberOfElts() / 2};
  il::Array<double> injection_rate{1, inj_rate};
  hfp2d::Sources Source(source_point, injection_rate);

  /////// INITIALIZATION OF SOLUTION OBJECT ///////
  // Solution object at time t_0, with pressure perturbation at t_0plus1
  il::Array<double> init_opening{2 * frac_mesh.numberOfElts(), 0.};
  il::Array<double> init_slip{2 * frac_mesh.numberOfElts(), 0.};
  const double err_frac_position = 2.;
  const double err_opening_dd = 2.;
  const double err_shear_dd = 2.;
  const double err_press = 2.;
  double shear_modulus = 1.;
  double shear_wave_vel = 1.;
  double damping_coeff = shear_modulus / (2 * shear_wave_vel);
  il::int_t init_iter_front_position = 0;
  il::int_t init_iter_ehls = 0;

  // Fault hydraulic diffusivity [L^2/T]
  double alpha = 10;
  il::Array<double> press_init_nodes{frac_mesh.numberOfNodes(), 0};
  // Pore pressure at time t_0plus1 at nodal points
  for (il::int_t j = 0; j < press_init_nodes.size(); ++j) {
    press_init_nodes[j] =
        inSituStress.getLocalInSituPorePressure(0) +
        (const_overpress *
         (erfc(fabs((frac_mesh.coordinates(j, 0)) / sqrt(alpha * t_ini)) / 2)));
  }

  // Pore pressure at time t_0plus1 at collocation points
  auto press_init_coll = il::dot(from_edge_to_coll_press, press_init_nodes);

  // Get the active set of collocation points at time t_0
  // by checking the Mohr-Coulomb criterion
  il::Array<int> init_failed_set_collpoints{0};
  init_failed_set_collpoints.reserve(2 * frac_mesh.numberOfElts());
  for (int j = 0, k = 0; j < 2 * frac_mesh.numberOfElts(); ++j) {
    if (inSituStress.uniformShearInSituTractions(frac_mesh)[j] >
        SolidEvolution.getFricCoeff(j) *
            (inSituStress.uniformNormalInSituTractions(frac_mesh)[j] -
             press_init_coll[j])) {
      init_failed_set_collpoints.resize(k + 1);
      init_failed_set_collpoints[k] = j;
      k = k + 1;
    }
  }

  // Get active set of elements at time t_0
  il::Array<int> init_set_elements{0};
  il::Array<int> outp{2};
  init_set_elements.reserve(2 * frac_mesh.numberOfElts());
  for (int l = 0, k = 0; l < init_failed_set_collpoints.size(); ++l, ++k) {
    init_set_elements.resize(k + 1);

    // Find 2D integer
    for (int i = 0; i < dof_single_dd.size(0); ++i) {
      for (int j = 0; j < dof_single_dd.size(1); ++j) {
        if (dof_single_dd(i, j) == init_failed_set_collpoints[l])
          outp[0] = i, outp[1] = j;
      }
    }

    init_set_elements[l] = outp[0];
  }

  // Delete duplicate elements
  il::Array<int> init_set_elmnts{};
  for (il::int_t i = 0; i < init_set_elements.size(); ++i) {
    bool already_there = false;
    for (il::int_t j = 0; j < init_set_elmnts.size(); ++j) {
      if (init_set_elements[i] == init_set_elmnts[j]) {
        already_there = true;
      }
    }
    if (!already_there) {
      init_set_elmnts.append(init_set_elements[i]);
    }
  }

  il::Array<int> init_active_set_elements{init_set_elmnts.size()};
  // Enforce the propagation to be element by element
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
      frac_mesh, t_ini, init_dt, init_opening, init_slip, press_init_nodes,
      inSituStress.uniformNormalInSituTractions(frac_mesh),
      inSituStress.uniformShearInSituTractions(frac_mesh),
      init_active_set_elements, init_iter_front_position, init_iter_ehls,
      err_frac_position, err_opening_dd, err_shear_dd, err_press);

  // Initialization of slippage length at time t_0
  // No slip as initial condition prior the pressurization -> Sl = 0
  double slipp_length_at_T0 = 0.;

  /////// LOOP IN TIME ///////
  std::string filename;
  double slipp_length_at_Tn = slipp_length_at_T0;
  double slipp_length_at_Tn_plus1;
  double current_crack_velocity;

  while ((SolutionAtTn.time() < t_max) &&
         (slipp_length_at_Tn <
                  abs(frac_mesh.coordinates(0, 0) -
           frac_mesh.coordinates(frac_mesh.numberOfNodes() - 1, 0))) &&
         (SolutionAtTn.frontIts() < SimulationParameters.frac_front_max_its)) {
    std::cout << " --------------------------------" << std::endl;
    std::cout << " Current time t = " << SolutionAtTn.time() << std::endl;

    hfp2d::Solution SolutionAtTnPlus1 =
        fractFrontPosition(kmat, from_edge_to_coll_dds, from_edge_to_coll_dd,
                           from_edge_to_coll_press, dof_single_dd, frac_mesh,
                           fracfluid, SimulationParameters, SolidEvolution,
                           FractureEvolution, Source, SolutionAtTn, expl_impl,
                           QD, damping_coeff, dilat_plast[0], inSituStress);

    // Calculate the new slippage length (i.e at time T_n+1)
    if (SolutionAtTnPlus1.activeElts().size() == 0) {
      slipp_length_at_Tn_plus1 = 0.;
    } else {
      slipp_length_at_Tn_plus1 =
          abs(frac_mesh.coordinates(SolutionAtTnPlus1.activeElts(0), 0) -
          frac_mesh.coordinates(SolutionAtTnPlus1.activeElts(
                                    SolutionAtTnPlus1.activeElts().size() - 1) +
                                    1,
                                0));
    }

    // Calculate the current crack velocity
    current_crack_velocity = (slipp_length_at_Tn_plus1 - slipp_length_at_Tn) /
                             SolutionAtTn.timestep();

    // Write solution at TnPlus1 to .Json file (solution which contains pressure
    // and active set of elements at TnPlus2 and so on..)
    filename = output_dir + std::string{"/Solution_at_time_"} +
               std::to_string(SolutionAtTn.time()) + std::string{".json"};
    SolutionAtTnPlus1.writeToFile(filename);

    // Set new time step and new time based on current crack velocity
    if (adapt_dt) {
      if (current_crack_velocity > 0.0) {
        double dt_new = 2 * frac_mesh.allEltSize()[0] / current_crack_velocity;
        if (dt_new > 3. * SolutionAtTnPlus1.timestep()) {
          SolutionAtTnPlus1.setTimeStep(3. * SolutionAtTnPlus1.timestep());
        } else {
          if (dt_new < 0.9 * SolutionAtTnPlus1.timestep()) {
            SolutionAtTnPlus1.setTimeStep(0.9 * SolutionAtTnPlus1.timestep());
          } else {
            SolutionAtTnPlus1.setTimeStep(dt_new);
          }
        }
      } else {
        SolutionAtTnPlus1.setTimeStep(dt_max);
      }
    }

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
    il::Array2D<il::int_t> &dof_single_dd, Mesh &theMesh,
    Fluid &FluidProperties, SimulationParameters &SimulationParameters,
    SolidEvolution &SolidEvolution, FractureEvolution &FractureEvolution,
    Sources &Source, Solution &SolutionAtTn, bool expl_impl, bool damping_term,
    double damping_coeff, double dilat_plast,
    InSituConditions &BackgroundLoadingConditions) {
  // Initialization of fracture front loop
  // SolutionAtTn -> solution object at current time Tn
  // SolutionAtTn_k -> current solution object at iter k of fracture front loop
  // SolutionAtTn_kPlus1 -> solution object at iter k+1 of fracture front loop
  hfp2d::Solution SolutionAtTn_k = SolutionAtTn;
  hfp2d::Solution SolutionAtTn_kPlus1;
  il::Status status;
  il::Norm norm;
  norm = il::Norm::L1;
  il::int_t cvg_front_posit = 0;   // convergence front position
  il::int_t iter_front_posit = 0;  // iteration fracture front position
  il::Array<int> dofs_act_elmts{};
  dofs_act_elmts.reserve(theMesh.numberDDDofs());
  il::Array<double> incrm_shearDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<double> incrm_openingDD{2 * theMesh.numberOfElts(), 0.};
  il::Array<int> new_active_set_elements{};
  new_active_set_elements.reserve(theMesh.numberOfElts());

  while ((cvg_front_posit != 1) &&
         (iter_front_posit < SimulationParameters.frac_front_max_its)) {
    ++iter_front_posit;

    std::cout << "\n  - Iter for fracture front position = " << iter_front_posit
              << std::endl;

    // Find the corresponding DOFs of the active elements
    for (il::int_t elmnt_i = 0, k = 0;
         elmnt_i < SolutionAtTn_k.activeElts().size(); ++elmnt_i) {
      dofs_act_elmts.resize(k + 4);
      for (il::int_t j = 0, l = k; j < theMesh.numberDDDofsPerElt(); ++j, ++l) {
        dofs_act_elmts[l] =
            (int)theMesh.dofDD(SolutionAtTn_k.activeElts(elmnt_i), j);
      }
      k = k + 4;
    }

    if (dofs_act_elmts.size() == 0) {
      std::cout
          << " SHEAR CRACK IS NOT ACTIVATED BUT PORE PRESSURE EVOLVE IN TIME !"
          << std::endl;
    }

    // Call the Reynolds solver for P1 elements
    SolutionAtTn_kPlus1 = hfp2d::reynoldsP1(
        theMesh, elast_matrix, fetc_dds, fetc_dd, fetc_press, SolutionAtTn,
        SolutionAtTn_k, incrm_shearDD, incrm_openingDD, SimulationParameters,
        FluidProperties, SolidEvolution, FractureEvolution, Source,
        dofs_act_elmts, status, norm, damping_term, damping_coeff, dilat_plast,
        BackgroundLoadingConditions);

    // Update active set of elements
    new_active_set_elements = SolutionAtTn_kPlus1.activeSetElements(
        theMesh, SolutionAtTn_kPlus1, SolidEvolution, fetc_press, dof_single_dd,
        SolutionAtTn_kPlus1.pressure());

    auto old_active_set_elements = SolutionAtTn_kPlus1.activeElts();

    if (new_active_set_elements.size() == 0) {
      cvg_front_posit = 1;
    } else {
      for (auto i : new_active_set_elements) {
        old_active_set_elements.append(i);
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
