//
// Created by lorenzo on 9/9/17.
//


#include "loadProperties.h"

namespace hfp2d {

Properties loadProperties(const Mesh &theLoadedMesh,
                          const il::String &inputFileName,
                          const il::MapArray<il::String, il::Dynamic> &propertiesMap) {
  //// Properties
  // We have to load properties for each material, and placing it into the respective slot in the Properties class.
  // For the moment we accept only one type of solid and fluid evolution.
  //
  // The loading works as following.
  // We start loading the elasticity properties which are -for the moment- constant on all the domain.
  double youngM = findDouble("Young_modulus",propertiesMap,inputFileName);
  double poissR = findDouble("Poisson_ratio",propertiesMap,inputFileName);
  ElasticProperties theSolid(youngM,poissR);

  // Then, we load the fluid properties assuming a newtonian fluid in the natural fractures.
  double flDens = findDouble("fluid_density",propertiesMap,inputFileName);
  double flCompr = findDouble("fluid_compressibility",propertiesMap,inputFileName);
  double flVisco = findDouble("fluid_viscosity",propertiesMap,inputFileName);
  Fluid theFluid(flDens,flVisco,flCompr);

  // Which type of solid evolution did we selected?
  il::String solidEvol = findString("solid_evolution_type", propertiesMap, inputFileName);
  il::String fluidEvol = findString("fluid_evolution_type", propertiesMap, inputFileName);

  /// ------------
  // HERE we should create a switch which (i) select the kind of constitutive model depending on the name
  // that is passed, (ii) allocates the memory for the parameters (iii) all stress/traction-separation laws/etc
  // are virtual functions which point to the correct ones.
  //
  // In addition, it should point to the correct loading procedure, since keywords can be different!
  //
  // The same should be done for the parameters and constitutive laws described by the fluid evolution types.
  //
  // For linear cohesive zone models we have failure stress and opening as parameters. Then we would like to save
  // also the last maximum opening. So we set it at 3.
  il::int_t solidEvolParamNum = 3;
  /// ------------


  // How many number of solid evolution are there?
  // We will maintain the same cohesive zone model for every element but with different parameters.
  il::int_t keyFound;
  il::int_t numMat = findInteger("number_of_materials",propertiesMap,inputFileName);

  if(theLoadedMesh.numMats() != numMat){
    std::cerr << "ERROR: mismatch between materials in geometry and materials in properties, \n"
              << "in file " << inputFileName << std::endl;
    exit(3);
  }

  //
  // Now we load the type of solid evolution.
  // For example: cohesive zone model with linear variation of the stresses w.r.t. aperture
  // The solid evolution is computed at each collocation point location, but the number of
  // historical variables and parameters are depending on the type of constitutive model.

  il::int_t numDisplDofs = theLoadedMesh.numDDDofs();
  il::int_t numPressDofs = theLoadedMesh.numPressDofs();

  il::Array<double> failureStresses(numDisplDofs);
  il::Array<double> maxOpenings(numDisplDofs);
  il::Array<double> permeabilities(numPressDofs);

  std::cout << "Num. Displ. DOFs " << numDisplDofs << std::endl;
  std::cout << "Num. Press. DOFs " << numPressDofs << std::endl;
  std::cout << "Num. Materials" << numMat << std::endl;
  std::cout << "Num. of Fractures " << theLoadedMesh.numFracs() << std::endl;
  std::cout << "Num. of Nodes " << theLoadedMesh.numNodes() << std::endl;
  std::cout << "Num. of Elements " << theLoadedMesh.numElems() << std::endl;
  std::cout << "Displ. DOF per Element " << theLoadedMesh.numDDDofsPerElem() << std::endl;
  std::cout << "Press. DOF per Element " << theLoadedMesh.numPressDofsPerElem()  << std::endl;
  std::cout << "Displ. DOF x Num. Elem." << theLoadedMesh.numDDDofsPerElem()* theLoadedMesh.numElems() << std::endl;
  std::cout << "Computed Num. of Displ. Dofs" << (theLoadedMesh.interpOrd()+1)* theLoadedMesh.numElems() << std::endl;
  std::cout << "Press. DOF x Num. Elem." << theLoadedMesh.numPressDofsPerElem()* theLoadedMesh.numElems() << std::endl;
  std::cout << "Computed Num. of Press. Dofs" << theLoadedMesh.numElems()* theLoadedMesh.interpOrd()+ theLoadedMesh.numFracs() << std::endl;

  // we scan along the vector of materials ID which nodes will have the material ID
  for(il::int_t materialID=0; materialID<numMat; materialID++){

    std::cout << "Material ID " << materialID << std::endl;

    // load the materialID-th material
    const il::String materialName = il::join("material", il::toString(materialID));

    keyFound = propertiesMap.search(materialName);

    if(propertiesMap.found(keyFound) && propertiesMap.value(keyFound).isMapArray()){

      // take the map array of the single material
      const il::MapArray<il::String, il::Dynamic> &singleMaterial = propertiesMap.value(keyFound).asMapArray();

      // load the parameters of the single material (in this case for the linear cohesive zone model)
      double singleFailureStress = findDouble("failure_stress", singleMaterial, inputFileName);
      double singleMaxOpening = findDouble("max_opening", singleMaterial, inputFileName);

      double singlePermeabiity = findDouble("permeability", singleMaterial, inputFileName);

      // save the loaded parameters only in those dofs which matID correspond to the one that has been loaded
      for(il::int_t elmtK=0; elmtK < theLoadedMesh.numElems(); elmtK++){
        if(theLoadedMesh.matID(elmtK)==materialID){

          std::cout << "Element " << elmtK << " mesh matID " << theLoadedMesh.matID(elmtK) << "materialID" << materialID << std::endl;

          // this loop is for collocation point properties (e.g. CZMs)
          for(il::int_t j=0; j< theLoadedMesh.numDDDofsPerElem(); j++){

            // save the material parameters at the location indicated by the dof handle
            failureStresses[theLoadedMesh.dofDD(elmtK, j)] = singleFailureStress;
            maxOpenings[theLoadedMesh.dofDD(elmtK, j)] = singleMaxOpening;

            //std::cout << j << " " << theLoadedMesh.dofDD(elmtK,j) << " " ;

          }
          //std::cout << std::endl;

          // this loop is for nodal properties (e.g. flow & transport)
          for(il::int_t j=0; j< theLoadedMesh.numPressDofsPerElem(); j++){

            permeabilities[theLoadedMesh.dofPress(elmtK,j)] = singlePermeabiity;
            //std::cout << j << " " << theLoadedMesh.dofPress(elmtK,j) << " " ;

          }
          //std::cout << std::endl;
        }
      }

    } else {

      std::cerr << "ERROR: missing Material number " << materialID << std::endl;
      std::cerr << "in file: " << inputFileName << std::endl;
      exit(4);

    }

  }
  std::cout << "Here " << failureStresses.size() << " " << maxOpenings.size() << " " << std::endl;

  SolidEvolution theSolidEvolution(failureStresses,maxOpenings);
  PermeabilityEvolution theFluidEvolution(permeabilities);

  std::cout << "Here " << std::endl;
  // Having loaded all the parameters, the final step is create the properties container from each single container

  Properties theProperties(theSolid,
                           theFluid,
                           theSolidEvolution,
                           theFluidEvolution);

  std::cout << "Here " << std::endl;
return theProperties;
}

}



