//
// Created by lorenzo on 9/14/17.
//

#include "SegmentData.h"


namespace hfp2d {

//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotation_matrix_2D( const double theta ) {
  il::StaticArray2D<double, 2, 2> R;

  R(0, 0) = cos(1. * theta);
  R(0, 1) = -1. * sin(1. * theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  return R;
}


}


/* ######################## OLD VERSION ###############################
//
// Created by lorenzo on 9/14/17.
//

#include "SegmentData.h"


namespace hfp2d {

//   Rotation Matrix
il::StaticArray2D<double, 2, 2> rotation_matrix_2D(double theta) {
  il::StaticArray2D<double, 2, 2> R;

  R(0, 0) = cos(1. * theta);
  R(0, 1) = -1. * sin(1. * theta);
  R(1, 0) = sin(theta);
  R(1, 1) = cos(theta);

  return R;
}

//

// Function returning the segment characteristic from the matrix of the
// coordinates of the end points and the knowledge of the degree of
// interpolation p
//  work for segment mesh
// Inputs
// mesh object
// ne element number in the mesh to get characterstic from
// p order of the interpolation for that mesh
SegmentData get_segment_DD_data(const Mesh &mesh,
                                il::int_t ne,
                                il::int_t p) {
  //  IL_ASSERT(Xs.size(0) == 2);
  //  IL_ASSERT(Xs.size(1) == 2);

  SegmentData segment;

  // compute element size
  il::StaticArray<double, 2> xdiff, s, n, xmean, xaux;
  il::Array2D<double> Xcol{p + 1, 2, 0};
  il::StaticArray2D<double, 2, 2> R;
  il::StaticArray2D<double, 2, 2> Xs;

  Xs(0, 0) = mesh.node(mesh.connectivity(ne, 0), 0);
  Xs(0, 1) = mesh.node(mesh.connectivity(ne, 0), 1);

  Xs(1, 0) = mesh.node(mesh.connectivity(ne, 1), 0);
  Xs(1, 1) = mesh.node(mesh.connectivity(ne, 1), 1);

  xdiff[0] = Xs(1, 0) - Xs(0, 0);
  xdiff[1] = Xs(1, 1) - Xs(0, 1);

  segment.size = sqrt(pow(xdiff[0], 2) + pow(xdiff[1], 2));

  // s=xdiff; // tangent vector
  s[0] = xdiff[0] / segment.size;
  s[1] = xdiff[1] / segment.size;
  n[0] = -1. * s[1];
  n[1] = s[0];  // normal vector

  segment.s = s;
  segment.n = n;

  segment.theta = acos(s[0] / sqrt(pow(s[0], 2) + pow(s[1], 2)));
  if (s[1] < 0) {
    segment.theta = -segment.theta;
  };

  // mid point of the element
  xmean[0] = (Xs(1, 0) + Xs(0, 0)) / 2.;
  xmean[1] = (Xs(1, 1) + Xs(0, 1)) / 2.;
  segment.Xmid = xmean;

  switch (p) {
  case 1: {  // linear DD
    Xcol(0, 0) = -1. / sqrt(2.);
    Xcol(0, 1) = 0.;
    Xcol(1, 0) = 1. / sqrt(2.);
    Xcol(1, 1) = 0.;
  };
    break;

  case 0: {
    Xcol(0, 0) = 0.;
    Xcol(0, 1) = 0.;
  };
    break;
  default:std::cout << "error\n";  //  error
    break;
  };

// Returning the collocation point in the global frame

  R = rotation_matrix_2D(segment.theta);

  for (int i = 0; i < p + 1; ++i) {
    xaux[0] = (segment.size) * Xcol(i, 0) / 2.;
    xaux[1] = (segment.size) * Xcol(i, 1) / 2.;

    xaux = il::dot(R, xaux);

    Xcol(i, 0) = xaux[0] + xmean[0];
    Xcol(i, 1) = xaux[1] + xmean[1];
  }

  segment.CollocationPoints = Xcol;

  return segment;  // return structure with all we need on the segment.
}
//----------------------------------------------------

}*/
