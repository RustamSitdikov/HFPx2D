//
// Created by lorenzo on 10/4/17.
//

#include "loadSources.h"

namespace hfp2d{

Sources loadSources(const Mesh &theLoadedMesh,
                    const il::String &inputFileName,
                    const il::MapArray<il::String, il::Dynamic> &sourcesMap) {

  il::int_t keyFound;

  // load number of states
  il::int_t numSources = findInteger("number_of_injection",sourcesMap,inputFileName);

  // size of array/vectors
  il::int_t numDisplDofs = theLoadedMesh.numDisplDofs();
  il::int_t numPressDofs = theLoadedMesh.numPressDofs();

  // creation of the array/vectors
  il::Array2D<double> injectionNodes(numSources,2);
  il::Array<double> injectionRates(numSources,0.0);

  // we scan along the vector of conditionID per element which element will have the condition ID
  for(il::int_t sourceID=0; sourceID<numSources; sourceID++){

    // load the sourceID-th value
    const il::String conditionName = il::join("source", il::toString(sourceID));

    keyFound = sourcesMap.search(conditionName);

    if(sourcesMap.found(keyFound) && sourcesMap.value(keyFound).isMapArray()){

      // take the map array of the single material
      const il::MapArray<il::String, il::Dynamic> &singleSources = sourcesMap.value(keyFound).asMapArray();

      // load the parameters of the single material (in this case for the linear cohesive zone model)
      double locationX = findDouble("location_X", singleSources, inputFileName);
      double locationY = findDouble("location_Y", singleSources, inputFileName);
      double injectionRate = findDouble("injection_rate", singleSources, inputFileName);

      injectionNodes(sourceID,0)=locationX;
      injectionNodes(sourceID,1)=locationY;

      injectionRates[sourceID]=injectionRate;

    } else {

      std::cerr << "ERROR: missing Condition number " << sourceID << std::endl;
      std::cerr << "in file: " << inputFileName << std::endl;
      exit(4);

    }

  }

  Sources theSources(injectionNodes,injectionRates);

  return theSources;


};

}