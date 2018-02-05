//
// Created by lorenzo on 9/18/17.
//


#include "LoadSingleMaterial.h"

namespace hfp2d {

/*
SolidEvolution loadSingleMaterial(const il::String &inputFileName,
                                  const il::int_t materialID,
                                  const il::MapArray<il::String, il::Dynamic> &propertiesMap) {

  il::String solidEvolType = findString("solid_evolution_type",propertiesMap,inputFileName);

  if (solidEvolType == "Cohesive Zone Model") {

    il::String czmType = findString("cohesive_zone_type",propertiesMap,inputFileName);

    if (czmType == "Pure Opening") {


    } else if (czmType == "Pure Shear") {

    } else if (czmType == "Mohr-Coulomb") {

    } else {
      std::cerr << "ERROR: unknown CZM type." << std::endl;
      std::cerr << "material: " << materialID << "file: " << inputFileName << std::endl;
      exit(3);
    }

  } else if (solidEvolType == "Continuum Damage") {

    // Placeholder for continuum damage

  } else if (solidEvolType == "Phase field") {

    // Placeholder for phase field modelling

  } else {

    std::cerr << "ERROR: unknown solid evolution type." << std::endl;
    std::cerr << "material: " << materialID << "file: " << inputFileName << std::endl;
    exit(3);

  }

  SolidEvolution theSolidEvol;

  return theSolidEvol;
}

*/

/*il::String findCZMType(const il::MapArray<il::String, il::Dynamic> &propertiesMap, const il::String &inputFileName) {

  il::int_t keyFound;
  il::String czmType;

  // load interpolation order of the mesh
  keyFound = propertiesMap.search("cohesive_zone_type");

  if (propertiesMap.found(keyFound) && propertiesMap.value(keyFound).isString()) {

    czmType = propertiesMap.value(keyFound).asString();

  } else {

    std::cerr << "ERROR: missing the CZM type in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return czmType;

}*/

/*il::String findSolidEvolutionType(const il::MapArray<il::String, il::Dynamic> &propertiesMap,
                                  const il::String &inputFileName) {

  il::int_t keyFound;
  il::String evolType;

  // load interpolation order of the mesh
  keyFound = propertiesMap.search("solid_evolution_type");

  if (propertiesMap.found(keyFound) && propertiesMap.value(keyFound).isString()) {

    evolType = propertiesMap.value(keyFound).asString();

  } else {

    std::cerr << "ERROR: missing the solid evolution in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return evolType;

}*/

}