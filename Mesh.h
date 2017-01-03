//
// Created by Brice Lecampion on 10.12.16.
//

#ifndef HFPX2D_MESH_H
#define HFPX2D_MESH_H

#include <il/Array2D.h>
#include <il/Array.h>
#include <il/StaticArray.h>
#include <il/StaticArray2D.h>


///// 1D mesh class
class Mesh {   // class for 1D mesh of 1D segment elements ?

 public:
  il::Array2D<double> Coor  ; // array of XY coordinates of nodes
  il::Array2D<int>  conn ;  // mesh connectivity array

  void set_values(il::Array2D<double>,il::Array2D<int>);

  int nelts() {return conn.size(0);} ;
  int ncoor() {return Coor.size(0);} ;

};


struct SegmentCharacteristic{
  double size;
  double theta; // angle w.r. to e_1
  il::StaticArray<double,2> n; // unit normal to segment in global system of coordinates
  il::StaticArray<double,2> s; // unit tangent to segment in global system of coordinates
  il::StaticArray<double,2> Xmid; //segment mid points coordinates.
  il::Array2D<double> CollocationPoints; //collocation points in global system of coordinates
};


il::StaticArray2D<double,2,2> rotation_matrix_2D(double theta);

SegmentCharacteristic get_segment_DD_characteristic(const il::Array2D<double> Xs , int const p);


#endif //HFPX2D_MESH_H
