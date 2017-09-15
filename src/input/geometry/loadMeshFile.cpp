//
// Created by lorenzo on 9/8/17.
//

#include "loadMeshFile.h"

namespace hfp2d {

void loadMeshFile(const il::String &meshFileName, il::io_t, Mesh &theMesh){

  // set a status variable
  il::Status status{};

  // load the file with a map array
  auto customMesh =
      il::load<il::MapArray<il::String, il::Dynamic>>(meshFileName, il::io, status);
  status.abortOnError();

  il::int_t keyFound;

  il::int_t numNodes;

  // Find the number of nodes
  keyFound = customMesh.search(il::toString("number_of_nodes"));
  if(customMesh.found(keyFound) && customMesh.value(keyFound).isInteger()){
   numNodes = customMesh.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: number of nodes not found." << std::endl;
    std::cerr << "In custom mesh file: " << meshFileName << std::endl;
    exit(2);
  }


  il::int_t numElements;

  // Find the number of elements
  keyFound = customMesh.search(il::toString("number_of_elements"));
  if(customMesh.found(keyFound) && customMesh.value(keyFound).isInteger()){
    numElements = customMesh.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: number of elements not found." << std::endl;
    std::cerr << "In custom mesh file: " << meshFileName << std::endl;
    exit(2);
  }


  il::Array2D<double> coordinates;

  // Find the coordinates
  keyFound = customMesh.search(il::toString("coordinates"));
  if(customMesh.found(keyFound) && customMesh.value(keyFound).isArray2dOfDouble()){
    coordinates = customMesh.value(keyFound).asArray2dOfDouble();
  } else {
    std::cerr << "ERROR: mesh coordinates not found." << std::endl;
    std::cerr << "In custom mesh file: " << meshFileName << std::endl;
    exit(2);
  }

  // Check consistency of data
  IL_EXPECT_FAST(coordinates.size(0) == numNodes);
  IL_EXPECT_FAST(coordinates.size(1) == 2);

/* TODO: part on hold due to missing implementation in IL
 * Already contacted the IL guys to add the missing definitions
 *
  il::Array2D<il::int_t> connectivity;

  // Find the coordinates
  keyFound = customMesh.search(il::toString("coordinates"));
  if(customMesh.found(keyFound) && customMesh.value(keyFound).isArray2dOfInt)){
    connectivity = customMesh.value(keyFound).asArray2dOfInt();
  } else {
    std::cerr << "ERROR: mesh connectivity not found." << std::endl;
    std::cerr << "In custom mesh file: " << meshFileName << std::endl;
    exit(2);
  }

  // Check consistency of data
  IL_EXPECT_FAST(connectivity.size(0) == numElements);
  IL_EXPECT_FAST(connectivity.size(1) == 2);


  il::Array<il::int_t> materialID;

  // Find the coordinates
  keyFound = customMesh.search(il::toString("material_ID"));
  if(customMesh.found(keyFound) && customMesh.value(keyFound).isArray2dOfInt()){
    materialID = customMesh.value(keyFound).asArrayOfInt();
  } else {
    std::cerr << "ERROR: material ID not found." << std::endl;
    std::cerr << "In custom mesh file: " << meshFileName << std::endl;
    exit(2);
  }

  // Check consistency of data
  IL_EXPECT_FAST(materialID.size() == numElements);


  il::Array<il::int_t> fractureID;

  // Find the coordinates
  keyFound = customMesh.search(il::toString("material_ID"));
  if(customMesh.found(keyFound) && customMesh.value(keyFound).isArray2dOfUint8()){
    fractureID = customMesh.value(keyFound).asArrayOfInt();
  } else {
    for(il::int i=0; i<numElements; i++){
      fractureID[i]=0;
    }
  }

  // Check consistency of data
  IL_EXPECT_FAST(fractureID.size() == numElements);

  // create the mesh object
  theMesh.init1DMesh(coordinates, connectivity, materialID, fractureID);

  */
}

}