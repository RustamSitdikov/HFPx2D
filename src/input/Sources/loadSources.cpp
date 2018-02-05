//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 05.02.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

#include "loadSources.h"

//
// Created by lorenzo on 10/4/17.
//

#include "loadSources.h"

namespace hfp2d {

Sources loadSources(const Mesh &theLoadedMesh, const il::String &inputFileName,
                    const il::MapArray<il::String, il::Dynamic> &sourcesMap) {
  il::int_t keyFound;

  // load number of states
  il::int_t numStates =
      findInteger("number_of_injection", sourcesMap, inputFileName);

  // size of array/vectors
  il::int_t numDisplDofs = theLoadedMesh.numberDDDofs();
  il::int_t numPressDofs = theLoadedMesh.numberPressDofs();

  // creation of the array/vectors
  il::Array<double> injectionRateVector(numPressDofs, 0.0);

  // we scan along the vector of conditionID per element which element will have
  // the condition ID
  for (il::int_t conditionID = 0; conditionID < numStates; conditionID++) {
    // load the conditionID-th value
    const il::String conditionName =
        il::join("source", il::toString(conditionID));

    keyFound = sourcesMap.search(conditionName);

    if (sourcesMap.found(keyFound) && sourcesMap.value(keyFound).isMapArray()) {
      // take the map array of the single material
      const il::MapArray<il::String, il::Dynamic> &singleSources =
          sourcesMap.value(keyFound).asMapArray();

      // load the parameters of the single material (in this case for the linear
      // cohesive zone model)
      double locationX = findDouble("location_X", singleSources, inputFileName);
      double locationY = findDouble("location_Y", singleSources, inputFileName);
      double injectionRate =
          findDouble("injection_rate", singleSources, inputFileName);

      il::int_t sourceNode =
          findSourceLocation(locationX, locationY, theLoadedMesh);

      injectionRateVector[sourceNode] = injectionRate;

    } else {
      std::cerr << "ERROR: missing Condition number " << conditionID
                << std::endl;
      std::cerr << "in file: " << inputFileName << std::endl;
      exit(4);
    }
  }

  Sources theSources(injectionRateVector);

  return theSources;
};

il::int_t findSourceLocation(const double locationX, const double locationY,
                             const Mesh &theLoadedMesh) {
  il::Array<double> squareOfLocation(theLoadedMesh.numberOfNodes());

  for (il::int_t i = 0; i < theLoadedMesh.numberOfNodes(); i++) {
    // double dist

    squareOfLocation[i] =
        sqrt(theLoadedMesh.node(i, 0) * theLoadedMesh.node(i, 0) +
             theLoadedMesh.node(i, 1) * theLoadedMesh.node(i, 1));
  }

  il::int_t indexOfMin = 0;
  for (il::int_t i = 1; i < theLoadedMesh.numberOfNodes(); i++) {
    // if(square)
  }

  return 0;
}
}