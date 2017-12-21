//
// HFPx2D project.
//
// Created by Brice Lecampion on 10.12.16.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//
//
//

// Inclusion from Inside Loop library
#include <il/Array2D.h>
#include <il/linear_algebra.h>
#include <il/math.h>

// Inclusion from the project
#include <src/core/ElasticProperties.h>
#include <src/core/Mesh.h>
#include "Simplified3D.h"

// contains fundamental plane-strain elasticity kernels.
// for fracture segment with linear variation of displacement discontinuities
// (DD)

namespace hfp2d {

//------------------------------------------------------------------------------
// implementation of the Wu & Olson (2015) Simplified 3D kernels for constant
// height
// fractures

// as per their notation, the heigth is along direction 2
// e_1 (x) corresponds to e_x in 2D plane-strain
// e_3 (z) corresponds to e_y in 2D plane-strain

// RONGVED SOLUTION FOR A P0 Rectangular dislocation in a full space
// dislocation is centered on the origin in the plane e_3=0 (z) , (-a,a) in e_1
// (x) (-b,b)
// in e_2 (z

// AUXILIARY FUNCTIONS NEEDED
// second order derivatives of I(x,y,z,xi,eta)
double ip11(double x, double y, double z, double xi, double eta) {
  //  (x - \[Xi])/((y - \[Eta] + Sqrt(Power(z,2) + Power(y - \[Eta],2) + Power(x
  //  - \[Xi],2)))*
  //      Sqrt(Power(z,2) + Power(y - \[Eta],2) + Power(x - \[Xi],2)))
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) / ((R + y - eta) * R);
}

double ip12(double x, double y, double z, double xi, double eta) {
  // double R ;

  return 1. / sqrt(x * x + y * y + z * z - 2 * y * eta + eta * eta -
                   2 * x * xi + xi * xi);
}

double ip13(double x, double y, double z, double xi, double eta) {
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + y - eta) * R);
}

double ip22(double x, double y, double z, double xi, double eta) {
  //  (y - \[Eta])/(Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (y - eta) / ((R + x - xi) * R);
}

double ip23(double x, double y, double z, double xi, double eta) {
  //  z/(Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + x - xi) * R);
}

double ip33(double x, double y, double z, double xi, double eta) {
  //  ((y - \[Eta]) (2 z^2 + (y - \[Eta])^2 + (x - \[Xi])^2) (x - \
//\[Xi]))/((z^2 + (y - \[Eta])^2) (z^2 + (x - \[Xi])^2) Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])
  double R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) * (y - eta) *
         (2 * z * z + (y - eta) * (y - eta) + (xi - x) * (xi - x)) /
         (R * (z * z + (x - xi) * (x - xi)) * (z * z + (y - eta) * (y - eta)));
}

//// third order derivatives of I(x,y,z,xi,eta)

double ip111(double x, double y, double z, double xi, double eta) {
  //  (R2 (Sqrt[R2] + y - \[Eta]) -
  //      Sqrt[R2] (x - \[Xi])^2 - (Sqrt[R2] + y - \[Eta]) (x - \[Xi])^2)/(R2^(
  //      3/2) (Sqrt[R2] + y - \[Eta])^2)

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * pow(x - xi, 2)) + (y - eta) * (R2 - pow(x - xi, 2))) /
         (pow(R2, 3. / 2.) * pow(R + y - eta, 2.));
}

double ip112(double x, double y, double z, double xi, double eta) {
  //  (-x + \[Xi])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (xi - x) / pow(R2, 3. / 2.);
}

double ip113(double x, double y, double z, double xi, double eta) {
  //-((z (y - \[Eta] +
  //  2 Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2]) (x - \[Xi]))/((y - \
//\[Eta] + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])^2 (z^2 + (y - \[Eta])^2 + \
//(x - \[Xi])^2)^(3/2)))

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return z * (xi - x) * (2. * R + y - eta) /
         (pow(R2, 3. / 2.) * pow(R + y - eta, 2));
}

double ip122(double x, double y, double z, double xi, double eta) {
  //(-y + \[Eta])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (eta - y) / pow(R2, 3. / 2.);
}

double ip123(double x, double y, double z, double xi, double eta) {
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return -z / pow(R2, 3. / 2.);
}

double ip133(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (y - \[Eta]))/(R2^(
  //      3/2) (R + y - \[Eta])^2)

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (R2 - z * z) * (y - eta)) /
         (pow(R2, 3. / 2.) * pow(R + y - eta, 2.));
}

double ip222(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 (y - \[Eta])^2) + (R2 - (y - \[Eta])^2) (x - \[Xi]))/(R2^(
  //      3/2) (R + x - \[Xi])^2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * pow(y - eta, 2.)) +
          (x - xi) * (R2 - (y - eta) * (y - eta))) /
         (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

double ip223(double x, double y, double z, double xi, double eta) {
  //  -((z (y - \[Eta]) (2 R + x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2))

  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return z * (eta - y) * (2 * R + x - xi) /
         (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

double ip233(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2)
  double R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  double R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (x - xi) * (R2 - z * z)) /
         (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

// CHEMMERY Integration function

typedef double (*vFunctionCall)(double x, double y, double z, double xi,
                                double eta);

double rectangular_integration(double x, double y, double z, double a, double b,
                               vFunctionCall Func) {
  double ma, mb;
  ma = -a;
  mb = -b;
  return (Func(x, y, z, a, b) - Func(x, y, z, a, mb) - Func(x, y, z, ma, b) +
          Func(x, y, z, ma, mb));
}

// Fundamental stress kernel - only  effect of DD_x(e_1) and DD_z (e3)
// this function output the 3D stress in the 3D frame for completeness ?
il::StaticArray2D<double, 2, 6> all_3d_stresses_kernel_s3d_p0_dd(
    double x, double y, double z, double a, double b, double G, double nu) {
  // x , y , z location where to compute stress
  //  a,b  1/2 size of the rectangular DD
  //  G Shear modulus, nu Poisson's ratio'
  // Ep instead of G ?
  //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
  //  DDx (shear), DDy (shear), DDz (normal)

  //  CONVENTION: Positive DDs in Overlap, Stresses positive in tension

  //  double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

  //  double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233;

  double Ce =
      1. * G / (4 * il::pi * (1. - nu));  // Minus sign here for convention !
  //  double sxx, sxy, sxz, syy, syz, szz;
  //
  il::StaticArray2D<double, 2, 6> Stress;
  // compute the Is function derivatives....

  double Ip11 = rectangular_integration(x, y, z, a, b, ip11);
  double Ip22 = rectangular_integration(x, y, z, a, b, ip22);
  double Ip33 = rectangular_integration(x, y, z, a, b, ip33);
  double Ip23 = rectangular_integration(x, y, z, a, b, ip23);
  double Ip12 = rectangular_integration(x, y, z, a, b, ip12);
  double Ip13 = rectangular_integration(x, y, z, a, b, ip13);

  double Ip111 = rectangular_integration(x, y, z, a, b, ip111);
  double Ip122 = rectangular_integration(x, y, z, a, b, ip122);
  double Ip133 = rectangular_integration(x, y, z, a, b, ip133);
  double Ip112 = rectangular_integration(x, y, z, a, b, ip112);
  double Ip113 = rectangular_integration(x, y, z, a, b, ip113);
  double Ip123 = rectangular_integration(x, y, z, a, b, ip123);
  //  Ip222 = rectangular_integration(x, y, z, a, b, ip222);
  double Ip233 = rectangular_integration(x, y, z, a, b, ip233);
  double Ip223 = rectangular_integration(x, y, z, a, b, ip223);

  // Stress row is dof (DDx,DDy,DDx), columns are sxx,syy,szz,sxy,sxz,syz

  // stress due to displacement discontinuity DDx (shear)
  Stress(0, 0) = Ce * (2. * Ip13 - z * Ip111);         // sxx
  Stress(0, 1) = Ce * (2. * nu * Ip13 - z * Ip122);    // syy
  Stress(0, 2) = Ce * (-z * Ip133);                    // szz
  Stress(0, 3) = Ce * ((1. - nu) * Ip23 - z * Ip112);  // sxy
  Stress(0, 4) = Ce * (Ip33 + nu * Ip22 - z * Ip113);  // sxz
  Stress(0, 5) = Ce * (-nu * Ip12 - z * Ip123);        // syz

  // stress du to displacement discontinuity  DDy (shear)
  //  Stress(1, 0) = Ce * (2. * nu * Ip23 - z * Ip112);    // sxx
  //  Stress(1, 1) = Ce * (2. * Ip23 - z * Ip222);         // syy
  //  Stress(1, 2) = Ce * (-z * Ip233);                    // szz
  //  Stress(1, 3) = Ce * ((1. - nu) * Ip13 - z * Ip122);  // sxy
  //  Stress(1, 4) = Ce * (-nu * Ip12 - z * Ip123);        // sxz
  //  Stress(1, 5) = Ce * (Ip33 + nu * Ip11 - z * Ip223);  // syz

  // stress du to displacement discontinuity DDz (normal)
  Stress(1, 0) = Ce * (Ip33 + (1. - 2. * nu) * Ip22 - z * Ip113);  // sxx
  Stress(1, 1) = Ce * (Ip33 + (1. - 2. * nu) * Ip11 - z * Ip223);  // syy
  Stress(1, 2) = Ce * (Ip33 - z * Ip113);                          // szz
  Stress(1, 3) = Ce * (-(1. - 2. * nu) * Ip12 - z * Ip123);        // sxy
  Stress(1, 4) = Ce * (-z * Ip133);                                // sxz
  Stress(1, 5) = Ce * (z * Ip233);                                 // syz

  return Stress;
}

//------------------------------------------------------------------------------

//  this function ouptuts the stress in the plane (x2=0) s11,s33,s13
//  simplified 3D kernel
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(double a, double b,
                                                          double G, double nu,
                                                          double xx,
                                                          double yy) {
  // xx,yy location where to compute stress in the 2D elastic plan
  //  a,b  1/2 size of the rectangular DD, b corresponds to the 1/2 FRACTURE
  //  HEIGHT
  //  G Shear modulus, nu Poisson's ratio'
  // Ep instead of G ?
  //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
  //  DDx (shear),  DDz (normal)
  //  CONVENTION: Positive DD in Overlap, stress positive in tension

  // switch to the 3D cartesian frame
  double x = xx;
  double z = yy;
  double y = 0.;

  double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

  double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233;

  double Ce = 1. * G / (4 * il::pi * (1. - nu));

  double C11, C12, D, H;

  D = sqrt(xx * xx + yy * yy);  // distance to DD center
  H = 2 * b;  // height -> may want to modify that & input H directly ?

  // correction factors
  C11 =
      1. - pow(D, 1.5) * pow(H / 1., 4.) / pow(D * D + H * H, (1.5 + 4.) / 2.);
  C12 =
      1. - pow(D, 1.8) * pow(H / 1., 4.) / pow(D * D + H * H, (1.8 + 4.) / 2.);

  //  double sxx, sxy, sxz, syy, syz, szz;
  //
  il::StaticArray2D<double, 2, 3> StressInPlane;
  // compute the Is function derivatives....

  //  Ip11 = rectangular_integration(x, y, z, a, b, ip11);
  Ip22 = rectangular_integration(x, y, z, a, b, ip22);
  Ip33 = rectangular_integration(x, y, z, a, b, ip33);
  //  Ip23 = rectangular_integration(x, y, z, a, b, ip23);
  //  Ip12 = rectangular_integration(x, y, z, a, b, ip12);
  Ip13 = rectangular_integration(x, y, z, a, b, ip13);

  Ip111 = rectangular_integration(x, y, z, a, b, ip111);
  //  Ip122 = rectangular_integration(x, y, z, a, b, ip122);
  Ip133 = rectangular_integration(x, y, z, a, b, ip133);
  //  Ip112 = rectangular_integration(x, y, z, a, b, ip112);
  Ip113 = rectangular_integration(x, y, z, a, b, ip113);
  //  Ip123 = rectangular_integration(x, y, z, a, b, ip123);
  //  Ip222 = rectangular_integration(x, y, z, a, b, ip222);
  //  Ip233 = rectangular_integration(x, y, z, a, b, ip233);
  //  Ip223 = rectangular_integration(x, y, z, a, b, ip223);

  // StressInPlane row is dof (DDshear,DDnormal),
  // columns are sxx,sxy,syy (in 2D plane frame)

  // stress due to displacement discontinuity DDx (shear)
  // sxx
  StressInPlane(0, 0) = Ce * C11 * (2. * Ip13 - z * Ip111);
  //  Stress(0, 1) = Ce * (2. * nu * Ip13 - z * Ip122);   // syy
  //  Stress(0, 3) = Ce * ((1. - nu) * Ip23 - z * Ip112);  // sxy
  // sxz (3) -> Correspond to sxy in 2d Plan elasticity frame
  StressInPlane(0, 1) = Ce * C12 * (Ip33 + nu * Ip22 - z * Ip113);
  //  Stress(0, 5) = Ce * (-nu * Ip12 - z * Ip123);        // syz
  // szz (3D)  -> correspond to syy in 2d Plan elasticity frame
  StressInPlane(0, 2) = Ce * C11 * (-z * Ip133);

  // stress du to displacement discontinuity DDz (normal)
  StressInPlane(1, 0) =
      Ce * C11 * (Ip33 + (1. - 2. * nu) * Ip22 - z * Ip113);  // sxx
  //  Stress(1, 3) = Ce * (-(1. - 2. * nu) * Ip12 - z * Ip123);        // sxy
  StressInPlane(1, 1) = Ce * C12 * (-z * Ip133);  // sxz
  //  Stress(1, 5) = Ce * (z * Ip233);                                 // syz
  //  Stress(1, 1) = Ce * (Ip33 + (1. - 2. * nu) * Ip11 - z * Ip223);  // syy
  StressInPlane(1, 2) = Ce * C11 * (Ip33 - z * Ip113);  // szz

  return StressInPlane;
}

//------------------------------------------------------------------------------

il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_s3d_dp0_dd(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t i_col,
    ElasticProperties Elas, double ker_options) {
  //  Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a piece-wise constant DD segment of size h  and height ker_options
  // (the   DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2])
  //
  // INPUTS:
  // St : returned stress
  // source_elt : element data structure of the source element
  // receiver_elt  : element data structure of the receiver element
  // i_col : integer for the collocation number where to compute the normal and
  // shear stress in the receiver element (0, 1)
  // Elas :: elastic properties object
  // ker_options : dummy argument here  for the height of the simplified 3D elt
  //
  // Note that we ensure .... that we return a 2*4 matrix, eventhough here
  // because
  // it's a piece-wise DD half is filled with zero

  // switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R =
      hfp2d::rotationMatrix2D(source_elt.theta());

  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = receiver_elt.CollocationPoints(i_col, i) - source_elt.Xmid(i);
  }
  xe = il::dot(R, xe);

  il::StaticArray<double, 2> n = il::dot(R, receiver_elt.n());
  il::StaticArray<double, 2> s = il::dot(R, receiver_elt.s());

  double h = source_elt.size();

  il::StaticArray2D<double, 2, 3> stress_l = stresses_kernel_s3d_p0_dd(
      h / 2., ker_options / 2., Elas.G(), Elas.nu(), xe[0], xe[1]);

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  il::StaticArray2D<double, 2, 4> St;

  // shear stress
  //  effect of shear dd
  St(0, 0) = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) +
             n2s2 * stress_l(0, 2);
  //  effect of normal dd
  St(0, 1) = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) +
             n2s2 * stress_l(1, 2);

  // normal stress
  //  effect of shear dd
  St(1, 0) = n1n1 * stress_l(0, 0) + 2. * n1n2 * stress_l(0, 1) +
             n2n2 * stress_l(0, 2);
  //  effect of normal dd
  St(1, 1) = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) +
             n2n2 * stress_l(1, 2);

  // padding zero due to this agnostic kernel assembly.... un-needed if working
  // by nodal values....

  St(0, 2) = 0;
  St(1, 2) = 0;
  St(0, 3) = 0;
  St(1, 3) = 0;

  return St;
};

//------------------------------------------------------------------------------
il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_s3d_dp0_dd_nodal(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t s_col,
    il::int_t i_col, ElasticProperties Elas, double ker_options) {
  //  Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a piece-wise constant DD segment of size h  and height ker_options
  // (the   DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2])
  //
  // INPUTS:
  // St : returned stress
  // source_elt : element data structure of the source element
  // receiver_elt  : element data structure of the receiver element
  // i_col : integer for the collocation number where to compute the normal and
  // shear stress in the receiver element (0, 1)
  // Elas :: elastic properties object
  // ker_options : dummy argument here  for the height of the simplified 3D elt
  //
  // Note that we ensure .... that we return a 2*4 matrix, eventhough here
  // because
  // it's a piece-wise DD half is filled with zero

  // switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R =
      hfp2d::rotationMatrix2D(source_elt.theta());

  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = receiver_elt.CollocationPoints(i_col, i) - source_elt.Xmid(i);
  }
  xe = il::dot(R, xe);

  il::StaticArray<double, 2> n = il::dot(R, receiver_elt.n());
  il::StaticArray<double, 2> s = il::dot(R, receiver_elt.s());

  double h = source_elt.size();

  il::StaticArray2D<double, 2, 3> stress_l = stresses_kernel_s3d_p0_dd(
      h / 2., ker_options / 2., Elas.G(), Elas.nu(), xe[0], xe[1]);

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  il::StaticArray2D<double, 2, 2> St;

  // shear stress
  //  effect of shear dd
  St(0, 0) = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) +
             n2s2 * stress_l(0, 2);
  //  effect of normal dd
  St(0, 1) = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) +
             n2s2 * stress_l(1, 2);

  // normal stress
  //  effect of shear dd
  St(1, 0) = n1n1 * stress_l(0, 0) + 2. * n1n2 * stress_l(0, 1) +
             n2n2 * stress_l(0, 2);
  //  effect of normal dd
  St(1, 1) = n1n1 * stress_l(1, 0) + 2. * n1n2 * stress_l(1, 1) +
             n2n2 * stress_l(1, 2);

  return St;
};
}
