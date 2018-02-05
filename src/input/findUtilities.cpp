//
// Created by lorenzo on 9/18/17.
//

#include "findUtilities.h"
namespace hfp2d {
il::String findString(const il::String &keyword,
                      const il::MapArray<il::String, il::Dynamic> &inputMap,
                      const il::String &inputFileName) {
  il::int_t foundKey;
  il::String foundString;

  // Look for keyword and store it in foundKey
  foundKey = inputMap.search(keyword);

  // First, check that the key is found
  if (inputMap.found(foundKey)) {

    if (inputMap.value(foundKey).isString()) {

      foundString = inputMap.value(foundKey).asString();

    } else {
      std::cerr << "ERROR: type mismatch for " << keyword
                << " in file " << inputFileName << std::endl;
      exit(3);
    }

  } else {

    std::cerr << "ERROR: missing keyword " << keyword
              << " in file " << inputFileName << std::endl;
    exit(3);

  }

  return foundString;

}

///////////////////////////////////////////////////////////////////////

il::int_t findInteger(const il::String &keyword,
                      const il::MapArray<il::String, il::Dynamic> &inputMap,
                      const il::String &inputFileName) {

  il::int_t foundKey;
  il::int_t foundInt;

  // Look for keyword and store it in foundKey
  foundKey = inputMap.search(keyword);

  // First, check that the key is found
  if (inputMap.found(foundKey)) {

    if (inputMap.value(foundKey).isInteger()) {

      foundInt = inputMap.value(foundKey).toInteger();

    } else {
      std::cerr << "ERROR: type mismatch for " << keyword
                << " in file " << inputFileName << std::endl;
      exit(3);
    }

  } else {

    std::cerr << "ERROR: missing keyword " << keyword
              << " in file " << inputFileName << std::endl;
    exit(3);

  }

  return foundInt;

}

///////////////////////////////////////////////////////////////////////

double findDouble(const il::String &keyword,
                  const il::MapArray<il::String, il::Dynamic> &inputMap,
                  const il::String &inputFileName) {

  il::int_t foundKey;
  double foundDouble;

  // Look for keyword and store it in foundKey
  foundKey = inputMap.search(keyword);

  // First, check that the key is found
  if (inputMap.found(foundKey)) {

    if (inputMap.value(foundKey).isDouble()) {

      foundDouble = inputMap.value(foundKey).toDouble();

    } else {
      std::cerr << "ERROR: type mismatch for " << keyword
                << " in file " << inputFileName << std::endl;
      exit(3);
    }

  } else {

    std::cerr << "ERROR: missing keyword " << keyword
              << " in file " << inputFileName << std::endl;
    exit(3);

  }

  return foundDouble;

};

///////////////////////////////////////////////////////////////////////

bool findBool(const il::String &keyword,
              const il::MapArray<il::String, il::Dynamic> &inputMap,
              const il::String &inputFileName) {

  il::int_t foundKey;
  bool foundBool;

  // Look for keyword and store it in foundKey
  foundKey = inputMap.search(keyword);

  // First, check that the key is found
  if (inputMap.found(foundKey)) {

    if (inputMap.value(foundKey).isBool()) {

      foundBool = inputMap.value(foundKey).toBool();

    } else {
      std::cerr << "ERROR: type mismatch for " << keyword
                << " in file " << inputFileName << std::endl;
      exit(3);
    }

  } else {

    std::cerr << "ERROR: missing keyword " << keyword
              << " in file " << inputFileName << std::endl;
    exit(3);

  }

  return foundBool;

}

///////////////////////////////////////////////////////////////////////

}