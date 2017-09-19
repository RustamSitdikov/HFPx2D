//
// Created by lorenzo on 9/9/17.
//


#include "loadProperties.h"

namespace hfp2d {

Properties loadProperties(const il::String &inputFileName,
                          const il::MapArray<il::String, il::Dynamic> &propertiesMap,
                          const il::int_t numOfMats) {

  Properties theProperties;
  //SolidEvolution
  LinearCZM theSolidEvolution;

  double youngModulus = findDouble("Young_modulus",propertiesMap,inputFileName);
  double poissonRatio = findDouble("Poisson_ratio",propertiesMap,inputFileName);

  il::String fluidLaw = findString("fluid_law",propertiesMap,inputFileName);
  double fluidDensity = findDouble("fluid_density",propertiesMap,inputFileName);
  double fluidCompres = findDouble("fluid_compressibility",propertiesMap,inputFileName);
  double fluidViscosity = findDouble("fluid_viscosity",propertiesMap,inputFileName);


  il::int_t keyFound;
  il::int_t numMaterials = findInteger("number_of_materials",propertiesMap,inputFileName);

  if(numOfMats != numMaterials) {
    std::cerr << "ERROR: mismatch in number of materials between mesh and properties." << std::endl;
    std::cerr << numOfMats << " in the mesh, " << numMaterials << "in the properties. file: " << inputFileName << std::endl;
    exit(2);
  }

  // For every fracture, search for the name "material" + number of material
  for (il::int_t materialID = 0; materialID < numMaterials; materialID++) {

    // now we create a string with the m name
    const il::String materialName = il::join("material", il::toString(materialID));

    // search for the material name
    keyFound = propertiesMap.search(materialName);

    // here we check that the "material#" is found and that it is a map array
    if (propertiesMap.found(keyFound) && propertiesMap.value(keyFound).isMapArray()) {

      // we save the map array into autoCreationMap
      const il::MapArray<il::String, il::Dynamic> &autoCreationMap = propertiesMap.value(keyFound).asMapArray();

      // and we pass it to the loader
      //theSolidEvolution = loadSingleMaterial(inputFileName, materialID, autoCreationMap);
      //theFluidEvolution = ...

    } else {
      std::cerr << "ERROR: missing material from list." << std::endl;
      std::cerr << "material: " << materialID << "file: " << inputFileName << std::endl;
      exit(2);
    }

  }


return theProperties;
}

/* old routines

double findYoungModulus(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName) {

  il::int_t keyFound;
  double youngModulus;

  // load interpolation order of the mesh
  keyFound = meshCreationMap.search("Young_modulus");

  if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isDouble()) {

    youngModulus = meshCreationMap.value(keyFound).toDouble();
    IL_EXPECT_FAST(youngModulus >= 0); // we do not want a negative Young's modulus

  } else {

    std::cerr << "ERROR: missing the Young modulus in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return youngModulus;

}

////////////////////////////////

double findPoissonRatio(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName) {

  il::int_t keyFound;
  double poissonRatio;

  // load interpolation order of the mesh
  keyFound = meshCreationMap.search("Poisson_ratio");

  if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isDouble()) {

    poissonRatio = meshCreationMap.value(keyFound).toDouble();
    IL_EXPECT_FAST(poissonRatio >= 0); // we do not want a negative Poisson ratio

  } else {

    std::cerr << "ERROR: missing the Poisson ratio in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return poissonRatio;

}

////////////////////////////////

double findFluidDensity(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName) {

  il::int_t keyFound;
  double fluidDensity;

  // load interpolation order of the mesh
  keyFound = meshCreationMap.search("fluid_density");

  if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isDouble()) {

    fluidDensity = meshCreationMap.value(keyFound).toDouble();
    IL_EXPECT_FAST(fluidDensity >= 0); // we do not want a negative fluid density

  } else {

    std::cerr << "ERROR: missing the fluid density in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return fluidDensity;

}

////////////////////////////////

double findFluidViscosity(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName) {

  il::int_t keyFound;
  double fluidViscosity;

  // load interpolation order of the mesh
  keyFound = meshCreationMap.search("fluid_viscosity");

  if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isDouble()) {

    fluidViscosity = meshCreationMap.value(keyFound).toDouble();
    IL_EXPECT_FAST(fluidViscosity >= 0); // we do not want a negative fluid viscosity

  } else {

    std::cerr << "ERROR: missing the fluid viscosity in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return fluidViscosity;

}

////////////////////////////////

double findFluidCompressibility(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName) {

  il::int_t keyFound;
  double fluidCompressibility;

  // load interpolation order of the mesh
  keyFound = meshCreationMap.search("fluid_compressibility");

  if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isDouble()) {

    fluidCompressibility = meshCreationMap.value(keyFound).toDouble();
    IL_EXPECT_FAST(fluidCompressibility >= 0); // we do not want a fluid compressibility

  } else {

    std::cerr << "ERROR: missing the fluid compressibility in properties." << std::endl;
    std::cerr << "file: " << inputFileName << std::endl;
    exit(3);

  }

  return fluidCompressibility;

}

  ////////////////////////////////

  il::int_t findNumMaterials(const il::MapArray<il::String, il::Dynamic> &meshCreationMap, const il::String &inputFileName) {

    il::int_t keyFound;
    il::int_t numMaterials;

    // load interpolation order of the mesh
    keyFound = meshCreationMap.search("number_of_materials");

    std::cout << meshCreationMap.found(keyFound) <<std::endl;
    std::cout << meshCreationMap.value(keyFound).isInteger() << std::endl;

    if (meshCreationMap.found(keyFound) && meshCreationMap.value(keyFound).isInteger()) {

      numMaterials = meshCreationMap.value(keyFound).toInteger();
      IL_EXPECT_FAST(numMaterials >= 0); // we do not want a negative number of materials

    } else {

      std::cerr << "ERROR: missing the number of materials in properties." << std::endl;
      std::cerr << "file: " << inputFileName << std::endl;
      exit(3);

    }

    return numMaterials;

  }
*/
}



