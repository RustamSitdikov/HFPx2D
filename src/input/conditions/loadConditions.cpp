//
// Created by lorenzo on 10/4/17.
//

#include "loadConditions.h"

namespace hfp2d {

Conditions loadConditions(const Mesh &theLoadedMesh,
                          const il::String &inputFileName,
                          const il::MapArray<il::String, il::Dynamic> &conditionsMap) {

  il::int_t keyFound;

  // load number of states
  il::int_t numStates = findInteger("number_of_conditions",conditionsMap,inputFileName);

  if(numStates=!theLoadedMesh.numberOfConditions()){
    std::cerr << "ERROR: mismatch between conditions in geometry and in-situ conditions, \n"
              << "in file " << inputFileName << std::endl;
    exit(4);
  }

  // size of array/vectors
  il::int_t numDisplDofs = theLoadedMesh.numberOfDisplDofs();
  il::int_t numPressDofs = theLoadedMesh.numberOfPressDofs();

  // creation of the array/vectors
  il::Array2D<double> stressDistribution(numDisplDofs,3);
  il::Array<double> porePressureDistribution(numPressDofs);

  // we scan along the vector of conditionID per element which element will have the condition ID
  for(il::int_t conditionID=0; conditionID<numStates; conditionID++){

    // load the conditionID-th value
    const il::String conditionName = il::join("condition", il::toString(conditionID));

    keyFound = conditionsMap.search(conditionName);

    if(conditionsMap.found(keyFound) && conditionsMap.value(keyFound).isMapArray()){

      // take the map array of the single material
      const il::MapArray<il::String, il::Dynamic> &singleCondition = conditionsMap.value(keyFound).asMapArray();

      // load the parameters of the single material (in this case for the linear cohesive zone model)
      double singleSxx = findDouble("s_xx", singleCondition, inputFileName);
      double singleSyy = findDouble("s_yy", singleCondition, inputFileName);
      double singleSxy = findDouble("s_xy", singleCondition, inputFileName);

      double singlePP = findDouble("pp", singleCondition, inputFileName);

      // save the loaded parameters only in those dofs which matID correspond to the one that has been loaded
      for(il::int_t elmtK=0; elmtK < numDisplDofs; elmtK++){
        if(theLoadedMesh.condID(elmtK)==conditionID){

          // this loop is for collocation point properties (e.g. CZMs)
          for(il::int_t j=0; j<theLoadedMesh.numberOfDisplDofsPerElement(); j++){

            // save the material parameters at the location indicated by the dof handle
            stressDistribution(theLoadedMesh.dofDispl(elmtK,j),0) = singleSxx;
            stressDistribution(theLoadedMesh.dofDispl(elmtK,j),1) = singleSyy;
            stressDistribution(theLoadedMesh.dofDispl(elmtK,j),2) = singleSxy;

          }

          // this loop is for nodal properties (e.g. flow & transport)
          for(il::int_t j=0; j<theLoadedMesh.numberOfPressDofsPerElement(); j++){

            porePressureDistribution[theLoadedMesh.dofPress(elmtK,j)] = singlePP;

          }

        }
      }

    } else {

      std::cerr << "ERROR: missing Condition number " << conditionID << std::endl;
      std::cerr << "in file: " << inputFileName << std::endl;
      exit(4);

    }

  }

  Conditions theConditions(stressDistribution,porePressureDistribution);

  return theConditions;

};

}