//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 29.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details.
//

// Inclusion from Inside Loop library
#include <il/Array2D.h>
#include <il/linear_algebra.h>
#include <il/math.h>

// Inclusion from the project
#include <src/core/ElasticProperties.h>
#include "PlaneStrainInfinite.h"

// contains fundamental plane-strain elasticity kernels.
// for fracture segment with linear variation of displacement discontinuities
// (DD)
// convention: Stress positive in tension
//      Displacement discontinuity positive in overlap dd=u^- - u^+
//  the equation to solve are schematically
//    n (sigma^o + sigma^induced ) n = t_n = - p (p fluid pressure, the minus is
//    because the normal is n=n^-)
//    s (sigma^o + sigma^induced ) n = t_s
//

namespace hfp2d {

il::StaticArray2D<double, 2, 4> stresses_kernel_dp1_dd(double h, double Ep,
                                                       double x, double y) {
  // function computing the stresses at (x,y) induced by a linear DD segment (of
  // total length h) centered on the origin [-h/2,h/2]
  // it returns stresses due to a linear variation from an unit value at the
  // left node (node 1)  to zero at the right node (coordinates 2) for both
  // shear and
  // opening displacement discontinuity
  // and stresses due to a linear variation from an unit value at the right
  // coordinates
  // (node 1) to zero at the right node (coordinates 2) for both shear and
  // opening
  // displacement discontinuity
  // Ep is the plane strain Young'smodulus
  // notation of stress component:
  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress
  // e.g.:  sxxs1 - sxx component due to a shear DD with a linear variation with
  // unit value at coordinates 1 .
  // note that we have the following relations : sxxn = sxys,  sxyn = syys such
  // that we don t output all values.

  //  should have a check for h>0
  double overhalfh = 2. / h;
  // minus sign here for convention of positive DD in overlap
  double Elascoef = -1. * Ep / (4. * il::pi);

  double xp = x * overhalfh;  // change of variable to perform the computation
  // for a unit segment [-1,1]
  double yp = y * overhalfh;  // change of variable to perform the computation
  // for a unit segment [-1,1]

  // auxiliary variables.
  double xpm1 = xp - 1.;
  double xpp1 = xp + 1.;
  double dxpp3 = 2. * xp + 3.;
  double xpm2 = xp - 2.;
  double r1 = xpm1*xpm1 + yp*yp;
  double r2 =  xpp1*xpp1 + yp*yp;
  double yp2m1 = yp*yp - 1;
  double yp2p1 = yp*yp + 1;

  double AtanAux1 = atanh((2. * xp) / (xp*xp  + yp2p1)) * 0.5;
  double AtanAux2 = atanh((2. * xp) / (xp*xp  + yp2p1)) * 0.5;

  double AtanAuxxpm1 = atan(xpm1 / yp);
  double AtanAuxxpp1 = atan(xpp1 / yp);

  double sxxs1, sxxs2;

  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress
// MY COMMMENT HERE 
  double sxys1 =
      AtanAux1 -
          (1./(r1*r2*r2)) * (
          (xpm1*xpm1) * pow(xpp1, 3) + 2. * xp * (3. + xp) * xpp1 * yp *yp + xpm1 * pow(yp, 4));

  double syys1 = (1./r1) *(1./(r2*r2))  * (2 * xpm1 * yp *  (xpp1*xpp1) -
                                                2 * (1 + 3 * xp) * pow(yp, 3));

  double syyn1 =
      AtanAux1 - xpp1 * (1/(r2*r2)) * ((xpp1*xpp1) + 3 * (yp*yp)) +
      2 * xp * (yp*yp) *
          pow(2 * yp2m1 * (xp*xp) + pow(xp, 4) + (yp2p1*yp2p1), -1);

  double sxys2 =
      -AtanAux2 +
      1/(r1*r1) * 1/(r2) *
          (pow(xpm1, 3) * pow(xpp1, 2) +
           2 * (-3 + xp) * xp * xpm1 * (yp*yp) + xpp1 * pow(yp, 4));

  double syys2 = pow(r1, -2) * pow(r2, -1) *
                 (2 * xpp1 * yp * pow(xpm1, 2) + (2 - 6 * xp) * pow(yp, 3));
  double syyn2 =
      -AtanAux2 +
      pow(r1, -2) * pow(r2, -1) *
          (pow(xpm1, 3) * pow(xpp1, 2) +
           2 * (2 + xp) * xpm1 * xpp1 * pow(yp, 2) + (-3 + xp) * pow(yp, 4));

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
  // row 1: effect of node 1, row 2 effect of coordinates 2
  // columns sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are
  // respectively equal to sxys and syys )

  il::StaticArray2D<double, 2, 4> Stress;

  // switch back to the segment [-h/2,h/2]
  //  we put directly here the constant -E'/(4Pi)
  Stress(0, 0) = Elascoef * overhalfh * sxxs1;
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

//------------------------------------------------------------------------------
// new api - for general kernel call for DDs
il::StaticArray2D<double, 2, 4> normal_shear_stress_kernel_dp1_dd(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t i_col,
    ElasticProperties Elas, double ker_options) {
  //   Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a linear DD segment of size h
  // (the linear DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2]
  // Material of Plane-strain Young's modulus Ep
  // INPUTS:
  // St : returned stress
  // source_elt : element data structure of the source element
  // receiver_elt  : element data structure of the receiver element
  // i_col : integer for the collocation point number where to compute the
  // normal and
  // shear stress in the receiver element (0, 1)
  // Elas :: elastic properties object
  // ker_options : dummy argument here (double) - needed for agnostic call ..

  // switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R =
      hfp2d::rotationMatrix2D(source_elt.theta());

  il::StaticArray2D<double, 2, 2> Rt=R;
  Rt(0,1)=R(1,0);
  Rt(1,0)=R(0,1);

  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = receiver_elt.CollocationPoints(i_col, i) - source_elt.Xmid(i);
  }

  xe = il::dot(Rt,xe);

  il::StaticArray<double, 2> n = il::dot(Rt, receiver_elt.n());
  il::StaticArray<double, 2> s = il::dot(Rt, receiver_elt.s());

  double h = source_elt.size();

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  // columns sxxs, sxys, syys, syyn
  // (knowing that sxxn and sxyn are respectively equal to sxys and syys )
  il::StaticArray2D<double, 2, 4> stress_l =
      stresses_kernel_dp1_dd(h, Elas.Ep(), xe[0], xe[1]);

  // shear stress
  // coordinates 1
  // shear dd
  double sh1s = n1s1 * stress_l(0, 0) + n1s2pn2s1 * stress_l(0, 1) +
                n2s2 * stress_l(0, 2);
  // normal dd
  double sh1n = n1s1 * stress_l(0, 1) + n1s2pn2s1 * stress_l(0, 2) +
                n2s2 * stress_l(0, 3);
  // coordinates 2
  // shear dd
  double sh2s = n1s1 * stress_l(1, 0) + n1s2pn2s1 * stress_l(1, 1) +
                n2s2 * stress_l(1, 2);
  // normal dd
  double sh2n = n1s1 * stress_l(1, 1) + n1s2pn2s1 * stress_l(1, 2) +
                n2s2 * stress_l(1, 3);

  // normal stress
  // coordinates 1
  // shear dd
  double sn1s =
      n1n1 * stress_l(0, 0) + 2 * n1n2 * stress_l(0, 1) + n2n2 * stress_l(0, 2);
  // normal dd
  double sn1n =
      n1n1 * stress_l(0, 1) + 2 * n1n2 * stress_l(0, 2) + n2n2 * stress_l(0, 3);
  // coordinates 2
  // shear dd
  double sn2s =
      n1n1 * stress_l(1, 0) + 2 * n1n2 * stress_l(1, 1) + n2n2 * stress_l(1, 2);
  // normal dd
  double sn2n =
      n1n1 * stress_l(1, 1) + 2 * n1n2 * stress_l(1, 2) + n2n2 * stress_l(1, 3);

  // output the desired stress components induced either by normal or shear dd
  // of the 2 nodes of the linear DD segment.
  il::StaticArray2D<double, 2, 4> St;
  // shear stress (node 1 shear dd, normal dd ; coordinates 2 shear dd , normal
  // dd)
  St(0, 0) = sh1s;
  St(0, 1) = sh1n;
  St(0, 2) = sh2s;
  St(0, 3) = sh2n;

  // normal stress (coordinates 1 shear dd, normal dd ; node 2 shear dd , normal
  // dd)
  St(1, 0) = sn1s;
  St(1, 1) = sn1n;
  St(1, 2) = sn2s;
  St(1, 3) = sn2n;
  // St dimensions 2 * 4

  return St;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//      KERNEL FUNCTION WORKING BY NODES - for hmat

il::StaticArray<double, 4> stresses_kernel_dp1_dd_nodal(il::int_t local_node_i,
                                                        double h, double Ep,
                                                        double x, double y) {
  // function computing the stresses at (x,y) induced by a linear DD segment (of
  // total length h) centered on the origin [-h/2,h/2]
  // it returns stresses due to a linear variation from an unit value at the
  // left node (node 1)  to zero at the right node (coordinates 2) for both
  // shear and
  // opening displacement discontinuity
  // OR stresses due to a linear variation from an unit value at the right
  // coordinates
  // (node 1) to zero at the right node (coordinates 2) for both shear and
  // opening
  // displacement discontinuity
  // Inputs
  // local_node_i :: integer describing if we need to compute the effect of the
  // left (1) or right (2) nodes
  // h :: element size
  // Ep ::  is the plane strain Young'smodulus
  // x :: x position of the point where to compute the stress (element frame)
  // y :: y position of the point where to compute the stress (element frame)

  // Notation of stress component:
  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress
  // e.g.:  sxxs1 - sxx component due to a shear DD with a linear variation with
  // unit value at coordinates 1 .
  // note that we have the following relations : sxxn = sxys,  sxyn = syys such
  // that we don t output all values.

  IL_EXPECT_FAST(local_node_i == 0 || local_node_i == 1);

  //  should have a check for h>0
  double overhalfh = 2. / h;
  // minus sign here for convention of positive DD in overlap
  double Elascoef = -1. * Ep / (4. * il::pi);

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

  double sxxs1, sxxs2;

  // 1 : left node influence, 2 : right coordinates influence
  // s : shear dof, n : normal dof
  // sxx xx stress, sxy xy stress, syy yy stress

  il::StaticArray<double, 4> Stress;

  if (local_node_i == 0) {  // node 1 (in C convention;))

    double sxys1 =
        AtanAux1 -
        pow(r1, -1.) * pow(r2, -2.) *
            (pow(xpm1, 2.) * pow(xpp1, 3.) +
             2. * xp * (3. + xp) * xpp1 * pow(yp, 2.) + xpm1 * pow(yp, 4.));

    double syys1 =
        pow(r1, -1.) * pow(r2, -2.) *
        (2 * xpm1 * yp * pow(xpp1, 2.) - 2 * (1 + 3 * xp) * pow(yp, 3.));

    double syyn1 =
        AtanAux1 - xpp1 * pow(r2, -2) * (pow(xpp1, 2) + 3 * pow(yp, 2)) +
        2 * xp * pow(yp, 2) *
            pow(2 * yp2m1 * pow(xp, 2) + pow(xp, 4) + pow(yp2p1, 2), -1);

    if (yp != 0.)  // case of observation point on the unit sgement line
    {
      sxxs1 = AtanAuxxpm1 - AtanAuxxpp1 +
              2 * yp * pow(r1, -1) * pow(r2, -2) *
                  (xpm1 * xpm2 * pow(xpp1, 2) + (3 + dxpp3 * xp) * pow(yp, 2) +
                   pow(yp, 4));  //

    } else {
      sxxs1 = ((il::pi) *
               (xpm1 / sqrt(pow(xpm1, 2.)) - xpp1 / sqrt(pow(xpp1, 2.)))) /
              2.;  // could be simplified further as 2 or 0 depending on abs(xp)
    }

    // switch back to the segment [-h/2,h/2]
    //  we put directly here the constant -E'/(4Pi)
    Stress[0] = Elascoef * overhalfh * sxxs1;  // Sxx due to shear opening
    Stress[1] = Elascoef * overhalfh * sxys1;  // Sxx due to normal opening
    Stress[2] = Elascoef * overhalfh * syys1;  // Syy due to shear opening
    Stress[3] = Elascoef * overhalfh * syyn1;
  };

  if (local_node_i == 1) {  // node 2 in C convention
    double sxys2 =
        -AtanAux2 +
        pow(r1, -2) * pow(r2, -1) *
            (pow(xpm1, 3) * pow(xpp1, 2) +
             2 * (-3 + xp) * xp * xpm1 * pow(yp, 2) + xpp1 * pow(yp, 4));

    double syys2 = pow(r1, -2) * pow(r2, -1) *
                   (2 * xpp1 * yp * pow(xpm1, 2) + (2 - 6 * xp) * pow(yp, 3));
    double syyn2 =
        -AtanAux2 +
        pow(r1, -2) * pow(r2, -1) *
            (pow(xpm1, 3) * pow(xpp1, 2) +
             2 * (2 + xp) * xpm1 * xpp1 * pow(yp, 2) + (-3 + xp) * pow(yp, 4));

    if (yp != 0.)  // case of observation point on the unit sgement line
    {
      sxxs2 = -AtanAuxxpm1 +
              pow(r1, -2) * pow(r2, -1) *
                  (AtanAuxxpp1 * r2 * pow(r1, 2) -
                   2 * (2 + xp) * xpp1 * yp * pow(xpm1, 2) -
                   2 * (3 + xp * (-3 + 2 * xp)) * pow(yp, 3) - 2 * pow(yp, 5));
    } else {
      sxxs2 = ((il::pi) *
               (-(xpm1 / sqrt(pow(xpm1, 2.))) + xpp1 / sqrt(pow(xpp1, 2.)))) /
              2.;
    }

    // switch back to the segment [-h/2,h/2]
    //  we put directly here the constant -E'/(4Pi)

    Stress[0] = Elascoef * overhalfh * sxxs2;
    Stress[1] = Elascoef * overhalfh * sxys2;
    Stress[2] = Elascoef * overhalfh * syys2;
    Stress[3] = Elascoef * overhalfh * syyn2;
  }

  // note that we have the following relations : sxxn = sxys,  sxyn = syys
  // we return a vector with 4 entries
  //   sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are
  // respectively equal to sxys and syys )

  return Stress;
}

////////////////////////////////////////////////////////////////////////////////
il::StaticArray2D<double, 2, 2> normal_shear_stress_kernel_dp1_dd_nodal(
    SegmentData source_elt, SegmentData receiver_elt, il::int_t s_col,
    il::int_t i_col, ElasticProperties Elas, double ker_options) {
  //   Function to get the normal and shear stress at a point on a surface
  // (with given normal and shear vector) induced by
  // a the local node s_col of a linear DD segment of size h
  // (the linear DD segment is assumed to be the unit element
  // along the x axis [-h/2,h/2]
  // Material of Plane-strain Young's modulus Ep
  // INPUTS:
  // St : returned stress
  // source_elt : element data structure of the source element
  // receiver_elt : element data structure of the receiver element
  // s_col : integer describing if we compute here the effect of node 0
  //        or 1 of the source element  (local node number here, i.e. 0 or 1)
  // i_col : integer for the collocation point number of the receiver element
  //        where to compute the normal and shear stress in the receiver element
  //        value either (0, 1)
  // Elas :: elastic properties object  (local node number here, i.e. 0 or 1)
  // ker_options : dummy argument here (double) - needed for agnostic call  due
  // to other kernel function
  // OUTPUT
  // 2*2 matrix containing the effect of node s_col of the source element, on
  // i_col
  // first row shear traction  (effect of shear and opening DD)
  // second row normal traction (effect of shear and opening DD)

  // switch to the frame of the source element....
  il::StaticArray2D<double, 2, 2> R =
      hfp2d::rotationMatrix2D(source_elt.theta());

  il::StaticArray2D<double, 2, 2> Rt=R;
  Rt(0,1)=R(1,0);
  Rt(1,0)=R(0,1);

  // receiver collocation point in the reference frame of the source element.
  il::StaticArray<double, 2> xe;
  for (int i = 0; i < 2; ++i) {
    xe[i] = receiver_elt.CollocationPoints(i_col, i) - source_elt.Xmid(i);
  }
  xe = il::dot(Rt, xe);

  il::StaticArray<double, 2> n = il::dot(Rt, receiver_elt.n());
  il::StaticArray<double, 2> s = il::dot(Rt, receiver_elt.s());

  double h = source_elt.size();

  double n1n1 = n[0] * n[0];
  double n2n2 = n[1] * n[1];
  double n1n2 = n[0] * n[1];
  double n1s1 = n[0] * s[0];
  double n2s2 = n[1] * s[1];

  double n1s2pn2s1 = n[0] * s[1] + n[1] * s[0];

  // columns sxxs, sxys, syys, syyn
  // (knowing that sxxn and sxyn are respectively equal to sxys and syys )
  il::StaticArray<double, 4> stress_l =
      stresses_kernel_dp1_dd_nodal(s_col, h, Elas.Ep(), xe[0], xe[1]);

  // shear stress
  // shear dd
  double shs =
      n1s1 * stress_l[0] + n1s2pn2s1 * stress_l[1] + n2s2 * stress_l[2];
  // normal dd
  double shn =
      n1s1 * stress_l[1] + n1s2pn2s1 * stress_l[2] + n2s2 * stress_l[3];

  // normal stress
  // shear dd
  double sns = n1n1 * stress_l[0] + 2 * n1n2 * stress_l[1] + n2n2 * stress_l[2];
  // normal dd
  double snn = n1n1 * stress_l[1] + 2 * n1n2 * stress_l[2] + n2n2 * stress_l[3];

  // output the desired stress components induced either by normal or shear dd
  // of the 2 nodes of the linear DD segment.
  il::StaticArray2D<double, 2, 2> St;

  // shear stress (node 1 shear dd, normal dd ; coordinates 2 shear dd , normal
  // dd)
  St(0, 0) = shs;
  St(0, 1) = shn;

  // normal stress (coordinates 1 shear dd, normal dd ; node 2 shear dd , normal
  // dd)
  St(1, 0) = sns;
  St(1, 1) = snn;

  return St;
}

//------------------------------------------------------------------------------
}