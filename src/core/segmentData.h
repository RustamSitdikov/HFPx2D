//
// Created by lorenzo on 9/14/17.
//

#ifndef HFPX2DUNITTEST_SEGMENTDATA_H
#define HFPX2DUNITTEST_SEGMENTDATA_H

#include <src/core/Mesh.h>
//#include <cmath>
//#include <il/Array.h>
//#include <il/Array2D.h>
//#include <il/StaticArray.h>
//#include <il/StaticArray2D.h>
//#include <il/container/1d/SmallArray.h>
//#include <il/String.h>
//#include <iostream>
#include <il/linear_algebra.h>



namespace hfp2d{

struct SegmentData {
  double size;
  double theta;  // angle w.r. to e_1
  // unit normal to segment in global system of coordinates
  il::StaticArray<double, 2> n;
  // unit tangent to segment in global system of coordinates
  il::StaticArray<double, 2> s;
  // segment mid points coordinates.
  il::StaticArray<double, 2> Xmid;
  // collocation points in global system of coordinates
  il::Array2D<double> CollocationPoints;
};

il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta);

SegmentData get_segment_DD_data(const Mesh &mesh,
                                il::int_t ne,
                                il::int_t p);


}

#endif //HFPX2DUNITTEST_SEGMENTDATA_H
