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

// contains fundamental plane-strain elasticity kernels.
// for fracture segment with linear variation of displacement discontinuities
// (DD)

#include <il/math.h>

#include "Elasticity2D.h"

namespace hfp2d {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(const double h,
                                                       const double Ep,
                                                       const double x,
                                                       const double y) {
  // function computing the stresses at (x,y) induced by a linear DD segment (of
  // total length h) centered on the origin [-h/2,h/2]
  // it returns stresses due to a linear variation from an unit value at the
  // left node (node 1)  to zero at the right node (node 2) for both shear and
  // opening displacement discontinuity
  // and stresses due to a linear variation from an unit value at the right node
  // (node 1) to zero at the right node (node 2) for both shear and opening
  // displacement discontinuit
  // Ep is the plane strain Young'smodulus
  // notation of stress component:
  // 1 : left node influence, 2 : right node influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress
  // e.g.:  sxxs1 - sxx component due to a shear DD with a linear variation with
  // unit value at node 1 .
  // note that we have the following relations : sxxn = sxys,  sxyn = syys such
  // that we don t output all values.

  //  should have a check for h>0
  double overhalfh = 2. / h;
  double Elascoef = Ep / (4. * il::pi);

  double xp = x * overhalfh;  // change of variable to perform the computation
                              // for a unit segment [-1,1]
  double yp = y * overhalfh;  // change of variable to perform the computation
                              // for a unit segment [-1,1]

  // auxiliary variables.
  double xpm1 = xp - 1.;
  double xpp1 = xp + 1.;
  double dxpp3 = 2. * xp + 3.;
  double xpm2 = xp - 2.;
  double r1 = pow(xpm1, 2.) + pow(yp, 2.);
  double r2 = pow(xpp1, 2.) + pow(yp, 2.);
  double yp2m1 = pow(yp, 2.) - 1;
  double yp2p1 = pow(yp, 2.) + 1;

  double AtanAux1 = atanh((2. * xp) / (pow(xp, 2.) + yp2p1)) * 0.5;
  double AtanAux2 = atanh((2. * xp) / (pow(xp, 2.) + yp2p1)) * 0.5;

  double AtanAuxxpm1 = atan(xpm1 / yp);
  double AtanAuxxpp1 = atan(xpp1 / yp);

  double sxxs1, sxys1, syys1, syyn1, sxxs2, sxys2, syyn2, syys2;

  // 1 : left node influence, 2 : right node influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress

  sxys1 = AtanAux1 -
          pow(r1, -1.) * pow(r2, -2.) *
              (pow(xpm1, 2.) * pow(xpp1, 3.) +
               2. * xp * (3. + xp) * xpp1 * pow(yp, 2.) + xpm1 * pow(yp, 4.));

  syys1 = pow(r1, -1.) * pow(r2, -2.) *
          (2 * xpm1 * yp * pow(xpp1, 2.) - 2 * (1 + 3 * xp) * pow(yp, 3.));

  syyn1 = AtanAux1 - xpp1 * pow(r2, -2) * (pow(xpp1, 2) + 3 * pow(yp, 2)) +
          2 * xp * pow(yp, 2) *
              pow(2 * yp2m1 * pow(xp, 2) + pow(xp, 4) + pow(yp2p1, 2), -1);

  sxys2 = -AtanAux2 +
          pow(r1, -2) * pow(r2, -1) *
              (pow(xpm1, 3) * pow(xpp1, 2) +
               2 * (-3 + xp) * xp * xpm1 * pow(yp, 2) + xpp1 * pow(yp, 4));

  syys2 = pow(r1, -2) * pow(r2, -1) *
          (2 * xpp1 * yp * pow(xpm1, 2) + (2 - 6 * xp) * pow(yp, 3));
  syyn2 = -AtanAux2 +
          pow(r1, -2) * pow(r2, -1) * (pow(xpm1, 3) * pow(xpp1, 2) +
                                       2 * (2 + xp) * xpm1 * xpp1 * pow(yp, 2) +
                                       (-3 + xp) * pow(yp, 4));

  if (yp != 0.)  // case of observation point on the unit sgement line
  {
    sxxs1 = AtanAuxxpm1 - AtanAuxxpp1 +
            2 * yp * pow(r1, -1) * pow(r2, -2) *
                (xpm1 * xpm2 * pow(xpp1, 2) + (3 + dxpp3 * xp) * pow(yp, 2) +
                 pow(yp, 4));  //

    sxxs2 = -AtanAuxxpm1 +
            pow(r1, -2) * pow(r2, -1) *
                (AtanAuxxpp1 * r2 * pow(r1, 2) -
                 2 * (2 + xp) * xpp1 * yp * pow(xpm1, 2) -
                 2 * (3 + xp * (-3 + 2 * xp)) * pow(yp, 3) - 2 * pow(yp, 5));
  } else {
    sxxs1 =
        ((il::pi) * (xpm1 / sqrt(pow(xpm1, 2.)) - xpp1 / sqrt(pow(xpp1, 2.)))) /
        2.;  // could be simplified further as 2 or 0 depending on abs(xp)
    sxxs2 = ((il::pi) *
             (-(xpm1 / sqrt(pow(xpm1, 2.))) + xpp1 / sqrt(pow(xpp1, 2.)))) /
            2.;
  }

  // note that we have the following relations : sxxn = sxys,  sxyn = syys
  // we return a matrix with 4 columns and 2 rows
  // row 1: effect of node 1, row 2 effect of node 2
  // columns sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are
  // respectively equal to sxys and syys )

  il::StaticArray2D<double, 2, 4> Stress;

  // switch back to the segment [-h/2,h/2]

  Stress(0, 0) = Elascoef * overhalfh *
                 sxxs1;  //  we put directly here the constant E'/(4Pi)
  Stress(0, 1) = Elascoef * overhalfh * sxys1;
  Stress(0, 2) = Elascoef * overhalfh * syys1;
  Stress(0, 3) = Elascoef * overhalfh * syyn1;

  Stress(1, 0) = Elascoef * overhalfh * sxxs2;
  Stress(1, 1) = Elascoef * overhalfh * sxys2;
  Stress(1, 2) = Elascoef * overhalfh * syys2;
  Stress(1, 3) = Elascoef * overhalfh * syyn2;

  return Stress;
}

//------------------------------------------------------------------------------
il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    const il::StaticArray<double, 2> xe, const double& h,
    const il::StaticArray<double, 2> s, const il::StaticArray<double, 2> n,
    const double Ep) {
  //   Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a linear DD segment of size h
  // (the linear DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2]
  // Material of Plane-strain Young's modulus Ep
  // INPUTS:
  // St:: returned stress
  // xe (x,y):: collocation point to compute the stress in the reference frame
  // of the unit element
  //
  //  h :: elt size
  // s :: tangent vector at xe on which to project to obtain the shear stress
  // n:: nornal vector at xe on which to project to obtain the normal stress

  double n1n1, n2n2, n1s1, n2s2, n1n2, n1s2pn2s1;

  il::StaticArray2D<double, 2, 4> stress_l;
  il::StaticArray2D<double, 2, 4> St;

  n1n1 = n[0] * n[0];
  n2n2 = n[1] * n[1];
  n1n2 = n[0] * n[1];
  n1s1 = n[0] * s[0];
  n2s2 = n[1] * s[1];

  n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  // columns sxxs, sxys, syys, syyn
  // (knowing that sxxn and sxyn are respectively equal to sxys and syys )
  stress_l = stresses_kernel_dp1_dd(h, Ep, xe[0], xe[1]);

  double sh1s, sn1s, sh1n, sn1n, sh2s, sn2s, sh2n, sn2n;
  // shear stress
  // node 1
  // shear dd
  sh1s = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) +
         n2s2 * stress_l(0, 2);
  // normal dd
  sh1n = n1s1 * stress_l(0, 1) + n1s2pn2s1 * stress_l(0, 2) +
         n2s2 * stress_l(0, 3);
  // node 2
  // shear dd
  sh2s = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) +
         n2s2 * stress_l(1, 2);
  // normal dd
  sh2n = n1s1 * stress_l(1, 1) + n1s2pn2s1 * stress_l(1, 2) +
         n2s2 * stress_l(1, 3);

  // normal stress
  // node 1
  // shear dd
  sn1s =
      n1n1 * stress_l(0, 0) + 2 * n1n2 * stress_l(0, 1) + n2n2 * stress_l(0, 2);
  // normal dd
  sn1n =
      n1n1 * stress_l(0, 1) + 2 * n1n2 * stress_l(0, 2) + n2n2 * stress_l(0, 3);
  // node 2
  // shear dd
  sn2s =
      n1n1 * stress_l(1, 0) + 2 * n1n2 * stress_l(1, 1) + n2n2 * stress_l(1, 2);
  // normal dd
  sn2n =
      n1n1 * stress_l(1, 1) + 2 * n1n2 * stress_l(1, 2) + n2n2 * stress_l(1, 3);

  // output the desired stress components induced either by normal or shear dd
  // of the 2 nodes of the linear DD segment.
  // shear stress (node 1 shear dd, normal dd ; node 2 shear dd , normal dd)
  St(0, 0) = sh1s;
  St(0, 1) = sh1n;
  St(0, 2) = sh2s;
  St(0, 3) = sh2n;
  // normal stress (node 1 shear dd, normal dd ; node 2 shear dd , normal dd)
  St(1, 0) = sn1s;
  St(1, 1) = sn1n;
  St(1, 2) = sn2s;
  St(1, 3) = sn2n;
  // St dimensions 2 * 4

  return St;
}

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
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) / ((R + y - eta) * R);
}

double ip12(double x, double y, double z, double xi, double eta) {
  // double R ;

  return 1. / sqrt(x * x + y * y + z * z - 2 * y * eta + eta * eta -
                   2 * x * xi + xi * xi);
}

double ip13(double x, double y, double z, double xi, double eta) {
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + y - eta) * R);
}

double ip22(double x, double y, double z, double xi, double eta) {
  //  (y - \[Eta])/(Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (y - eta) / ((R + x - xi) * R);
}

double ip23(double x, double y, double z, double xi, double eta) {
  //  z/(Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] (x + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2] - \[Xi]))
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return z / ((R + x - xi) * R);
}

double ip33(double x, double y, double z, double xi, double eta) {
  //  ((y - \[Eta]) (2 z^2 + (y - \[Eta])^2 + (x - \[Xi])^2) (x - \
//\[Xi]))/((z^2 + (y - \[Eta])^2) (z^2 + (x - \[Xi])^2) Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])
  double R;
  R = sqrt((x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z);

  return (x - xi) * (y - eta) *
         (2 * z * z + (y - eta) * (y - eta) + (xi - x) * (xi - x)) /
         (R * (z * z + (x - xi) * (x - xi)) * (z * z + (y - eta) * (y - eta)));
}

//// third order derivatives of I(x,y,z,xi,eta)

double ip111(double x, double y, double z, double xi, double eta) {
  //  (R2 (Sqrt[R2] + y - \[Eta]) -
  //      Sqrt[R2] (x - \[Xi])^2 - (Sqrt[R2] + y - \[Eta]) (x - \[Xi])^2)/(R2^(
  //      3/2) (Sqrt[R2] + y - \[Eta])^2)

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * pow(x - xi, 2)) + (y - eta) * (R2 - pow(x - xi, 2))) /
         (pow(R2, 3. / 2.) * pow(R + y - eta, 2.));
}

double ip112(double x, double y, double z, double xi, double eta) {
  //  (-x + \[Xi])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (xi - x) / pow(R2, 3. / 2.);
}

double ip113(double x, double y, double z, double xi, double eta) {
  //-((z (y - \[Eta] +
  //  2 Sqrt[z^2 + (y - \[Eta])^2 + (x - \[Xi])^2]) (x - \[Xi]))/((y - \
//\[Eta] + Sqrt[
  //  z^2 + (y - \[Eta])^2 + (x - \[Xi])^2])^2 (z^2 + (y - \[Eta])^2 + \
//(x - \[Xi])^2)^(3/2)))

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return z * (xi - x) * (2. * R + y - eta) /
         (pow(R2, 3. / 2.) * pow(R + y - eta, 2));
}

double ip122(double x, double y, double z, double xi, double eta) {
  //(-y + \[Eta])/(z^2 + (y - \[Eta])^2 + (x - \[Xi])^2)^(3/2)
  double R2;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return (eta - y) / pow(R2, 3. / 2.);
}

double ip123(double x, double y, double z, double xi, double eta) {
  double R2;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;

  return -z / pow(R2, 3. / 2.);
}

double ip133(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (y - \[Eta]))/(R2^(
  //      3/2) (R + y - \[Eta])^2)

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * z * z) + (R2 - z * z) * (y - eta)) /
         (pow(R2, 3. / 2.) * pow(R + y - eta, 2.));
}

double ip222(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 (y - \[Eta])^2) + (R2 - (y - \[Eta])^2) (x - \[Xi]))/(R2^(
  //      3/2) (R + x - \[Xi])^2)
  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return (R * (R2 - 2. * pow(y - eta, 2.)) +
          (x - xi) * (R2 - (y - eta) * (y - eta))) /
         (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

double ip223(double x, double y, double z, double xi, double eta) {
  //  -((z (y - \[Eta]) (2 R + x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2))

  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

  return z * (eta - y) * (2 * R + x - xi) /
         (pow(R2, 3. / 2.) * pow(R + x - xi, 2.));
}

double ip233(double x, double y, double z, double xi, double eta) {
  //  (R (R2 - 2 z^2) + (R2 - z^2) (x - \[Xi]))/(R2^(3/2) (R + x - \[Xi])^2)
  double R2, R;
  R2 = (x - xi) * (x - xi) + (y - eta) * (y - eta) + z * z;
  R = sqrt(R2);

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
    double& x, double& y, double& z, double& a, double& b, double& G,
    double& nu) {
  // x , y , z location where to compute stress
  //  a,b  1/2 size of the rectangular DD
  //  G Shear modulus, nu Poisson's ratio'
  // Ep instead of G ?
  //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
  //  DDx (shear), DDy (shear), DDz (normal)

  double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

  double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233;

  double Ce = G / (4 * il::pi * (1. - nu));
  //  double sxx, sxy, sxz, syy, syz, szz;
  //
  il::StaticArray2D<double, 2, 6> Stress;
  // compute the Is function derivatives....

  Ip11 = rectangular_integration(x, y, z, a, b, ip11);
  Ip22 = rectangular_integration(x, y, z, a, b, ip22);
  Ip33 = rectangular_integration(x, y, z, a, b, ip33);
  Ip23 = rectangular_integration(x, y, z, a, b, ip23);
  Ip12 = rectangular_integration(x, y, z, a, b, ip12);
  Ip13 = rectangular_integration(x, y, z, a, b, ip13);

  Ip111 = rectangular_integration(x, y, z, a, b, ip111);
  Ip122 = rectangular_integration(x, y, z, a, b, ip122);
  Ip133 = rectangular_integration(x, y, z, a, b, ip133);
  Ip112 = rectangular_integration(x, y, z, a, b, ip112);
  Ip113 = rectangular_integration(x, y, z, a, b, ip113);
  Ip123 = rectangular_integration(x, y, z, a, b, ip123);
  //  Ip222 = rectangular_integration(x, y, z, a, b, ip222);
  Ip233 = rectangular_integration(x, y, z, a, b, ip233);
  Ip223 = rectangular_integration(x, y, z, a, b, ip223);

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

//  this function ouptuts the stress in the plane (x2=0) s11,s33,s13
//  simplified 3D kernel
il::StaticArray2D<double, 2, 3> stresses_kernel_s3d_p0_dd(
    const double a, const double b, const double G, const double nu,
    const double xx, const double yy) {
  // xx,yy location where to compute stress in the 2D elastic plan
  //  a,b  1/2 size of the rectangular DD, b corresponds to the 1/2 FRACTURE
  //  HEIGHT
  //  G Shear modulus, nu Poisson's ratio'
  // Ep instead of G ?
  //  Rectangular DDon plan z=0   x \in [-a,a], y \in [-b,b]
  //  DDx (shear),  DDz (normal)

  // switch to the 3D cartesian frame
  double x = xx;
  double z = yy;
  double y = 0.;

  double Ip11, Ip22, Ip33, Ip23, Ip12, Ip13;

  double Ip111, Ip122, Ip133, Ip112, Ip113, Ip123, Ip222, Ip223, Ip233;

  double Ce = G / (4 * il::pi * (1. - nu));

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

il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_s3d_dp0_dd(
    const il::StaticArray<double, 2>& xe, const double hx, const double height,
    const il::StaticArray<double, 2>& s, const il::StaticArray<double, 2>& n,
    const double G, const double nu) {
  il::StaticArray2D<double, 2, 3> stress_l;

  stress_l =
      stresses_kernel_s3d_p0_dd(hx / 2., height / 2., G, nu, xe[0], xe[1]);

  double n1n1, n2n2, n1s1, n2s2, n1n2, n1s2pn2s1;

  il::StaticArray2D<double, 2, 2> St;

  n1n1 = n[0] * n[0];
  n2n2 = n[1] * n[1];
  n1n2 = n[0] * n[1];
  n1s1 = n[0] * s[0];
  n2s2 = n[1] * s[1];

  n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

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
