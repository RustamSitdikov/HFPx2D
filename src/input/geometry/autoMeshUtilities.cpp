//
// Created by lorenzo on 9/13/17.
//

#include "autoMeshUtilities.h"


// Find X of center of mesh
double findXC(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap){

  il::int_t keyFound;

  double x_c;
  keyFound = autoCreationMap.search("x_c");
  if (autoCreationMap.found(keyFound)) {
    x_c = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_c in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return x_c;
}


// Find Y of center of mesh
double findYC(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap){

  il::int_t keyFound;

  double y_c;
  keyFound = autoCreationMap.search("y_c");
  if (autoCreationMap.found(keyFound)) {
    y_c = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_c in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return y_c;

}


double findX1(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  double x_1;
  keyFound = autoCreationMap.search(il::toString("x_1"));
  if (autoCreationMap.found(keyFound)) {
    x_1 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return x_1;
}

double findY1(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;
  double y_1;
  keyFound = autoCreationMap.search(il::toString("y_1"));
  if (autoCreationMap.found(keyFound)) {
    y_1 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return y_1;
}


double findX2(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;
  double x_2;
  keyFound = autoCreationMap.search(il::toString("x_2"));
  if (autoCreationMap.found(keyFound)) {
    x_2 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing x_1 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return x_2;
}

double findY2(const il::String &inputFileName,
              const il::int_t fractureID,
              const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;
  double y_2;
  keyFound = autoCreationMap.search(il::toString("y_2"));
  if (autoCreationMap.found(keyFound)) {
    y_2 = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing y_2 in custom automatic mesh." << std::endl;
    std::cerr << "layer:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return y_2;
}



// Find angle of diagonal mesh
double findAngle(const il::String &inputFileName,
                 il::int_t fractureID,
                 const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  double angle;
  keyFound = autoCreationMap.search("angle");
  if (autoCreationMap.found(keyFound)) {
    angle = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing angle in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return angle;
}


// Find length of mesh
double findLength(const il::String &inputFileName,
                  const il::int_t fractureID,
                  const il::MapArray<il::String, il::Dynamic> &autoCreationMap){

  il::int_t keyFound;

  double length;
  keyFound = autoCreationMap.search("length");
  if (autoCreationMap.found(keyFound)) {
    length = autoCreationMap.value(keyFound).toDouble();
  } else {
    std::cerr << "ERROR: missing length in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return length;

}


// Find number of elements for the mesh
il::int_t findNumElem(const il::String &inputFileName,
                      const il::int_t fractureID,
                      const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  il::int_t numElements;
  keyFound = autoCreationMap.search("number_of_elements");
  if (autoCreationMap.found(keyFound)) {
    numElements = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the number of elements in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return numElements;
}


// Find the interpolation order for the mesh
il::int_t findInterpOrder(const il::String &inputFileName,
                          const il::int_t fractureID,
                          const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  il::int_t numElements;
  keyFound = autoCreationMap.search("interpolation_order");
  if (autoCreationMap.found(keyFound)) {
    numElements = autoCreationMap.value(keyFound).toInteger();
  } else {
    std::cerr << "ERROR: missing the interpolation order in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);
  }

  return numElements;
}


// Find the location of the injection
il::String findSource(const il::String &inputFileName,
                      il::int_t fractureID,
                      const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  // active source of injection
  il::String sourceLocation;
  keyFound = autoCreationMap.search("active_injection");

  // if active injection keyword is found
  if (autoCreationMap.found(keyFound)) {

    // and active injection value is true
    if(autoCreationMap.value(keyFound).toBool()){ // active injection

      // look for the keyword location
      keyFound = autoCreationMap.search("location");

      // if the keyword location is found
      if(autoCreationMap.found(keyFound)){

        // then save the location of the active source (which is either start, center, end)
        sourceLocation = autoCreationMap.value(keyFound).asString();

      } else { // if location is not found

        // default location is the center BUT we make it mandatory to write it
        std::cerr << "ERROR: source location is missing in automatic mesh." << std::endl;
        std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
        exit(2);

      }

    } else { // if the active injection value is false

      // no active injection in this fault
      sourceLocation = "None";

    }

  } else { // if source term is missing

    std::cerr << "ERROR: source missing in automatic mesh." << std::endl;
    std::cerr << "fracture:" << fractureID << ", file: " << inputFileName << std::endl;
    exit(2);

  }

  return sourceLocation;
}


// Find the material ID of the mesh
il::int_t findMaterialID(const il::String &inputFileName,
                         const il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  // Material ID is not necessary. If it is not set, then it is zero.
  il::int_t materialID;
  // if the material_ID keyword is found
  keyFound = autoCreationMap.search("material_ID");
  if (autoCreationMap.found(keyFound)) {
    // save its value
    materialID = autoCreationMap.value(keyFound).toInteger();
  } else {
    // save zero as default value
    materialID = 0;
  }

  return materialID;
}


// Find the far field stress ID for the mesh
il::int_t findFarFieldID(const il::String &inputFileName,
                         const il::int_t fractureID,
                         const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  // farFieldStress ID is looked for, but if not found is default
  il::int_t farFieldID;
  keyFound = autoCreationMap.search("far_field_stress_ID");
  if (autoCreationMap.found(keyFound)) {
    // if the far_field_stress_ID keyword is found
    farFieldID = autoCreationMap.value(keyFound).toInteger();
  } else {
    // save zero as default value
    farFieldID = 0;
  }

  return farFieldID;
}


// Find the pore pressure ID for the mesh
il::int_t findPorePresID(const il::String &inputFileName,
                        const il::int_t fractureID,
                        const il::MapArray<il::String, il::Dynamic> &autoCreationMap) {

  il::int_t keyFound;

  // porePressCond ID is taken into consideration as farFieldStress
  il::int_t porePresID;
  // if the pore_pressure_ID keyword is found
  keyFound = autoCreationMap.search("pore_pressure_ID");
  if (autoCreationMap.found(keyFound)) {
    porePresID = autoCreationMap.value(keyFound).toInteger();
  } else {
    // save zero as default value
    porePresID = 0;
  }

  return porePresID;
}

