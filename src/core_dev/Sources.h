//
// Created by lorenzo on 10/4/17.
//

#ifndef HFPX2D_SOURCES_H
#define HFPX2D_SOURCES_H

#include <il/Array.h>


// todo : will evolve .... we may not need a class for that.
namespace hfp2d {
class Sources {

private:
  il::Array<il::int_t> injection_nodes_;
  il::Array<double> injection_rates_;

public:
  Sources() {};

  Sources(il::Array<il::int_t> injNodes,
          il::Array<double> injRates) {

    // we need the corresponding location of the injection
    // element number (for P0) or nodes number for P1 ....
    //

    // can be a list if there is multiple rate....

    this->injection_nodes_ = injNodes;
    this->injection_rates_ = injRates;


  };

  Sources(il::Array2D<double> injLocation,
          il::Array<double> injRates){



  }

};


il::int_t findSourceLocation(const double locationX,
                             const double locationY,
                             const Mesh &theLoadedMesh){

  il::Array<double> squareOfLocation(theLoadedMesh.numNodes());

  for(il::int_t i=0; i< theLoadedMesh.numNodes(); i++){

    double distanceOnX = locationX - theLoadedMesh.node(i,0);
    double distanceOnY = locationY - theLoadedMesh.node(i,1);


    squareOfLocation[i]=sqrt(distanceOnX*distanceOnX+
        distanceOnY*distanceOnY);

  }

  auto smallestDistance = std::max_element(std::begin(squareOfLocation),std::end(squareOfLocation));

  il::int_t indexOfMin=std::distance(std::begin(squareOfLocation), smallestDistance);

/*
  for(il::int_t i=1; i<theLoadedMesh.numNodes(); i++){
    if(squareOfLocation[i]<squareOfLocation[indexOfMin]){
      indexOfMin=i;
    }
  }
*/

  return indexOfMin;
}


}
#endif //HFPX2D_SOURCES_H
