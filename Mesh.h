//
// Created by Brice Lecampion on 10.12.16.
//

#ifndef HFPX2D_MESH_H
#define HFPX2D_MESH_H


#include <il/Array2D.h>
///// 1D mesh class
class Mesh {   // class for 1D mesh of 1D segment elements ?

 public:
  il::Array2D<double> Coor  ; // array of XY coordinates of nodes
  il::Array2D<int>  conn ;  // mesh connectivity array

  void set_values(il::Array2D<double>,il::Array2D<int>);

  int nelts() {return conn.size(0);} ;
  int ncoor() {return Coor.size(0);} ;

};

#endif //HFPX2D_MESH_H
