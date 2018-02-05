//
// Created by lorenzo on 9/9/17.
//

#include "LoadProperties.h"

namespace hfp2d {

void loadProperties(const Mesh &theLoadedMesh, const il::String &input_filename,
                    const il::MapArray<il::String, il::Dynamic> &propertiesMap,
                    ElasticProperties &ElasticProperties,
                    FluidProperties &FluidProperties,
                    SolidEvolution &SolidEvolution,
                    FractureEvolution &FractureEvolution) {
  
  //// Elastic Properties
  // We have to load the elastic properties for each material, and placing it
  // into the respective slot in the Properties class.
  // For the moment we accept only one type of solid and fluid evolution.
  //
  // The loading works as following.
  // We start loading the elasticity Properties which are -for the moment-
  // constant on all the domain.
  double young_modulus =
      findDouble("young_modulus", propertiesMap, input_filename);
  double poisson_ratio =
      findDouble("poisson_ratio", propertiesMap, input_filename);

  ElasticProperties.setElastProp(young_modulus, poisson_ratio);

  // Then, we load the fluid Properties assuming a newtonian fluid in the
  // natural fractures.
  double fluidDensity =
      findDouble("fluid_density", propertiesMap, input_filename);
  double fluidCompres =
      findDouble("fluid_compressibility", propertiesMap, input_filename);
  double fluidViscosity =
      findDouble("fluid_viscosity", propertiesMap, input_filename);

  FluidProperties.setFluidProp(fluidDensity, fluidViscosity, fluidCompres);

  // How many number of solid evolution are there?
  // We will maintain the same cohesive zone model for every element but with
  // different parameters.
  il::int_t keyFound;
  il::int_t numMaterials =
      findInteger("number_of_materials", propertiesMap, input_filename);

  // TODO: check that number_of_material = size of vector of material_ID that we
  // input in the .toml file

  // Now we load the type of solid evolution & fracture evolution.

  // Initialization
  il::Array<double> vector_peak_fric_coeff(theLoadedMesh.numberDDDofs() / 2);
  il::Array<double> vector_resid_fric_coeff(theLoadedMesh.numberDDDofs() / 2);
  il::Array<double> vector_resid_slip(theLoadedMesh.numberDDDofs() / 2);

  il::Array<double> vector_singleFailureStress(theLoadedMesh.numberDDDofs() /
                                               2);
  il::Array<double> vector_singleMaxOpening(theLoadedMesh.numberDDDofs() / 2);
  il::Array<double> vector_singlePermeability(theLoadedMesh.numberDDDofs() /
                                              2);

  il::Array<double> vector_initial_permeab(theLoadedMesh.numberDDDofs() / 2);
  il::Array<double> vector_increment_permeab(theLoadedMesh.numberDDDofs() /
                                             2);
  il::Array<double> vector_initial_hydr_width(theLoadedMesh.numberDDDofs() /
                                              2);
  il::Array<double> vector_increment_hydr_width(theLoadedMesh.numberDDDofs() /
                                                2);

  il::Array2D<il::int_t> dof_single_dd{theLoadedMesh.numberOfElts(),
                                       (theLoadedMesh.interpolationOrder() + 1),
                                       0};
  for (il::int_t i = 0; i < theLoadedMesh.numberOfElts(); i++) {
    for (il::int_t j = 0; j < 1 * (theLoadedMesh.interpolationOrder() + 1);
         j++) {
      dof_single_dd(i, j) =
          i * 1 * (theLoadedMesh.interpolationOrder() + 1) + j;
    }
  }

  // We scan along the vector of materials ID which nodes will have the material
  // ID
  for (il::int_t materialID = 1; materialID <= numMaterials; materialID++) {
    // load the materialID-th material
    const il::String materialName =
        il::join("material", il::toString(materialID));

    keyFound = propertiesMap.search(materialName);

    if (propertiesMap.found(keyFound) &&
        propertiesMap.value(keyFound).isMapArray()) {
      // take the map array of the single material
      const il::MapArray<il::String, il::Dynamic> &singleMaterial =
          propertiesMap.value(keyFound).asMapArray();

      // load the parameters of the single material (in this case for the linear
      // cohesive zone model)

      /// SOLID EVOLUTION TYPE
      il::String solid_evol_type =
          findString("solid_evolution_type", singleMaterial, input_filename);

      if (solid_evol_type == "Cohesive Zone Model") {
        il::String cohesive_zone_type =
            findString("cohesive_zone_type", singleMaterial, input_filename);

        if (cohesive_zone_type == "Pure Shear") {
          double peak_fric_coeff = findDouble("peak_friction_coefficient",
                                              singleMaterial, input_filename);

          double resid_fric_coeff = findDouble("residual_friction_coefficient",
                                               singleMaterial, input_filename);

          double resid_slip =
              findDouble("residual_slip", singleMaterial, input_filename);

          for (il::int_t elmt_k = 0; elmt_k < theLoadedMesh.numberOfElts();
               ++elmt_k) {
            if (theLoadedMesh.matid(elmt_k) == materialID) {
              for (il::int_t j = 0; j < theLoadedMesh.numberDDDofsPerElt() / 2;
                   ++j) {
                vector_peak_fric_coeff[dof_single_dd(elmt_k, j)] =
                    peak_fric_coeff;

                vector_resid_fric_coeff[dof_single_dd(elmt_k, j)] =
                    resid_fric_coeff;
                vector_resid_slip[dof_single_dd(elmt_k, j)] = resid_slip;
              }
            }
          }

        } else if (cohesive_zone_type == "Pure Opening") {
          // Put whathever you need for 'pure opening'
          double singleFailureStress =
              findDouble("failure_stress", singleMaterial, input_filename);
          double singleMaxOpening =
              findDouble("max_opening", singleMaterial, input_filename);
          double singlePermeability =
              findDouble("permeability", singleMaterial, input_filename);

          for (il::int_t elmt_k = 0; elmt_k < theLoadedMesh.numberOfElts();
               ++elmt_k) {
            if (theLoadedMesh.matid(elmt_k) == materialID) {
              for (il::int_t j = 0; j < theLoadedMesh.numberDDDofsPerElt() / 2;
                   ++j) {
                vector_singleFailureStress[dof_single_dd(elmt_k, j)] =
                    singleFailureStress;
                vector_singleMaxOpening[dof_single_dd(elmt_k, j)] =
                    singleMaxOpening;
                vector_singlePermeability[dof_single_dd(elmt_k, j)] =
                    singlePermeability;
              }
            }
          }

        } else if (cohesive_zone_type == "Opening and Shear") {
          double peak_fric_coeff = findDouble("peak_friction_coefficient",
                                              singleMaterial, input_filename);

          double resid_fric_coeff = findDouble("residual_friction_coefficient",
                                               singleMaterial, input_filename);

          double resid_slip =
              findDouble("residual_slip", singleMaterial, input_filename);

          double singleFailureStress =
              findDouble("failure_stress", singleMaterial, input_filename);
          double singleMaxOpening =
              findDouble("max_opening", singleMaterial, input_filename);
          double singlePermeability =
              findDouble("permeability", singleMaterial, input_filename);

          for (il::int_t elmt_k = 0; elmt_k < theLoadedMesh.numberOfElts();
               ++elmt_k) {
            if (theLoadedMesh.matid(elmt_k) == materialID) {
              for (il::int_t j = 0; j < theLoadedMesh.numberDDDofsPerElt() / 2;
                   ++j) {
                vector_peak_fric_coeff[dof_single_dd(elmt_k, j)] =
                    peak_fric_coeff;
                vector_resid_fric_coeff[dof_single_dd(elmt_k, j)] =
                    resid_fric_coeff;
                vector_resid_slip[theLoadedMesh.dofDD(elmt_k, j)] = resid_slip;
                vector_singleFailureStress[dof_single_dd(elmt_k, j)] =
                    singleFailureStress;
                vector_singleMaxOpening[dof_single_dd(elmt_k, j)] =
                    singleMaxOpening;
                vector_singlePermeability[dof_single_dd(elmt_k, j)] =
                    singlePermeability;
              }
            }
          }

        } else {
          std::cerr << "ERROR: wrong input for cohesive_zone_type in material "
                    << materialID << std::endl;
          std::cerr << "in file: " << input_filename << std::endl;
          std::cerr << "Possible inputs for material " << materialID
                    << " must be: 'Pure Opening', 'Pure Shear' or "
                       "'Opening and Shear' "
                    << materialID << std::endl;
          exit(EXIT_FAILURE);
        }
      } else {
        std::cerr << "ERROR: wrong solid_evolution_type in material "
                  << materialID << std::endl;
        std::cerr << "in file: " << input_filename << std::endl;
        std::cerr << "The input for material " << materialID
                  << " must be: 'Cohesive Zone Model' " << std::endl;
        exit(EXIT_FAILURE);
      }

      /// FRACTURE EVOLUTION TYPE
      il::String fracture_evol_type =
          findString("fracture_evolution_type", singleMaterial, input_filename);

      if (fracture_evol_type == "Pure Shear") {
        double initial_permeab =
            findDouble("initial_permeability", singleMaterial, input_filename);
        double increment_permeab = findDouble("increment_permeability",
                                              singleMaterial, input_filename);
        double initial_hydr_width = findDouble("initial_hydraulic_width",
                                               singleMaterial, input_filename);
        double increment_hydr_width = findDouble(
            "increment_hydraulic_width", singleMaterial, input_filename);

        for (il::int_t elmt_k = 0; elmt_k < theLoadedMesh.numberOfElts();
             ++elmt_k) {
          if (theLoadedMesh.matid(elmt_k) == materialID) {
            for (il::int_t j = 0; j < theLoadedMesh.numberDDDofsPerElt() / 2;
                 ++j) {
              vector_initial_permeab[dof_single_dd(elmt_k, j)] =
                  initial_permeab;
              vector_increment_permeab[dof_single_dd(elmt_k, j)] =
                  increment_permeab;
              vector_initial_hydr_width[dof_single_dd(elmt_k, j)] =
                  initial_hydr_width;

              vector_increment_hydr_width[dof_single_dd(elmt_k, j)] =
                  increment_hydr_width;
            }
          }
        }

      } else if (fracture_evol_type == "Pure opening") {
        double initial_permeab =
            findDouble("initial_permeability", singleMaterial, input_filename);
        double initial_hydr_width = findDouble("intitial_hydraulic_width",
                                               singleMaterial, input_filename);

        for (il::int_t elmt_k = 0; elmt_k < theLoadedMesh.numberOfElts();
             ++elmt_k) {
          if (theLoadedMesh.matid(elmt_k) == materialID) {
            for (il::int_t j = 0; j < theLoadedMesh.numberDDDofsPerElt() / 2;
                 ++j) {
              vector_initial_permeab[dof_single_dd(elmt_k, j)] =
                  initial_permeab;
              vector_initial_hydr_width[dof_single_dd(elmt_k, j)] =
                  initial_hydr_width;
            }
          }
        }

      } else {
        std::cerr << "ERROR: wrong fracture_evolution_type in material "
                  << materialID << std::endl;
        std::cerr << "in file: " << input_filename << std::endl;
        std::cerr << "Possible inputs for material " << materialID
                  << " must be: 'Pure Shear' or 'Pure Opening' " << std::endl;
        exit(EXIT_FAILURE);
      }

    } else {
      std::cerr << "ERROR: missing Material number " << materialID << std::endl;
      std::cerr << "in file: " << input_filename << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // Call the constructors for SolidEvolution and FractureEvolution
  SolidEvolution.setSolidEvolution(
      vector_peak_fric_coeff, vector_peak_fric_coeff, vector_resid_fric_coeff,
      vector_resid_slip, vector_singleFailureStress, vector_singleMaxOpening);

  FractureEvolution.setFractureEvolution(
      vector_initial_permeab, vector_increment_permeab, vector_resid_slip,
      vector_initial_hydr_width, vector_increment_hydr_width);
}
}
