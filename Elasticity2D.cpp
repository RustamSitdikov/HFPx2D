//
// Created by Brice Lecampion on 10.12.16.
//
// contains fundamental plane-strain elasticity kernels.
// for fracture segment with linear variation of displacement discontinuities (DD)

#include "Elasticity2D.h"

#include <il/math.h>


il::StaticArray2D<double,2,4> StressesKernelLinearDD(const double h,const double Ep,const double x,const double y)
{
// function computing the stresses at (x,y) induced by a linear DD segment (of total length h) centered on the origin [-h/2,h/2]
// it returns stresses due to a linear variation from an unit value at the left node (node 1)  to zero at the right node (node 2) for both shear and opening displacement discontinuity
// and stresses due to a linear variation from an unit value at the right node (node 1) to zero at the right node (node 2) for both shear and opening displacement discontinuit
// Ep is the plane strain Young'smodulus
// notation of stress component:
// 1 : left node influence, 2 : right node influence
// s : shear dof, n : normal dof
// sxx xx stress, sxy xy stress, syy yy stress
// e.g.:  sxxs1 - sxx component due to a shear DD with a linear variation with unit value at node 1 .
// note that we have the following relations : sxxn = sxys,  sxyn = syys such that we don t output all values.

//  should have a check for h>0
  double overhalfh =  2./h;
  double Elascoef = Ep / (4.*il::pi);

  double xp = x*overhalfh ;   // change of variable to perform the computation for a unit segment [-1,1]
  double yp = y*overhalfh;    // change of variable to perform the computation for a unit segment [-1,1]

// auxiliary variables.
  double xpm1 = xp -1.;
  double xpp1 = xp+1.;
  double dxpp3 =2.*xp+3.;
  double xpm2= xp-2.;
  double r1 = pow(xpm1,2.)+pow(yp,2.);
  double r2 = pow(xpp1,2.)+pow(yp,2.);
  double yp2m1=pow(yp,2.)-1;
  double yp2p1=pow(yp,2.)+1;

  double AtanAux1= atanh((2.*xp)/(pow(xp,2.) + yp2p1)) * 0.5 ;
  double AtanAux2= atanh((2.*xp)/(pow(xp,2.) + yp2p1)) * 0.5 ;

  double AtanAuxxpm1 =atan(xpm1/yp);
  double AtanAuxxpp1 =atan(xpp1/yp);

  double sxxs1,sxys1,syys1,syyn1,sxxs2,sxys2,syyn2,syys2 ;

// 1 : left node influence, 2 : right node influence
// s : shear dof, n : normal dof
// sxx xx stress, sxy xy stress, syy yy stress

  sxys1 =AtanAux1 - pow(r1,-1.)*pow(r2,-2.)*(pow(xpm1,2.)*pow(xpp1,3.) +
      2.*xp*(3.+ xp)*xpp1*pow(yp,2.) + xpm1*pow(yp,4.)) ;

  syys1 =pow(r1,-1.)*pow(r2,-2.)*(2*xpm1*yp*pow(xpp1,2.) - 2*(1 + 3*xp)*pow(yp,3.)) ;

  syyn1 = AtanAux1 - xpp1*pow(r2,-2)*(pow(xpp1,2) + 3*pow(yp,2)) +
      2*xp*pow(yp,2)*pow(2*yp2m1*pow(xp,2) + pow(xp,4) + pow(yp2p1,2),-1) ;

  sxys2 = -AtanAux2 + pow(r1,-2)*pow(r2,-1)*(pow(xpm1,3)*pow(xpp1,2) +
      2*(-3 + xp)*xp*xpm1*pow(yp,2) + xpp1*pow(yp,4));

  syys2 = pow(r1,-2)*pow(r2,-1)*(2*xpp1*yp*pow(xpm1,2) + (2 - 6*xp)*pow(yp,3));
  syyn2=-AtanAux2 + pow(r1,-2)*pow(r2,-1)*(pow(xpm1,3)*pow(xpp1,2) +
      2*(2 + xp)*xpm1*xpp1*pow(yp,2) + (-3 + xp)*pow(yp,4)) ;

  if (yp!=0.) // case of observation point on the unit sgement line
  {

    sxxs1 = AtanAuxxpm1 - AtanAuxxpp1 + 2*yp*pow(r1,-1)*pow(r2,-2)*(xpm1*xpm2*pow(xpp1,2) + (3 + dxpp3*xp)*pow(yp,2) + pow(yp,4)) ;//

    sxxs2=-AtanAuxxpm1 + pow(r1,-2)*pow(r2,-1)*
        (AtanAuxxpp1*r2*pow(r1,2) - 2*(2 + xp)*xpp1*yp*pow(xpm1,2) -
            2*(3 + xp*(-3 + 2*xp))*pow(yp,3) - 2*pow(yp,5));
  }
  else
  {
    sxxs1 = ((il::pi)*(xpm1/sqrt(pow(xpm1,2.)) - xpp1/sqrt(pow(xpp1,2.))))/2.;    //could be simplified further as 2 or 0 depending on abs(xp)
    sxxs2 = ((il::pi)*(-(xpm1/sqrt(pow(xpm1,2.))) + xpp1/sqrt(pow(xpp1,2.))))/2.;
  }

// note that we have the following relations : sxxn = sxys,  sxyn = syys
// we return a matrix with 4 columns and 2 rows
// row 1: effect of node 1, row 2 effect of node 2
// columns sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are respectively equal to sxys and syys )

  il::StaticArray2D<double,2,4> Stress  ;

  // switch back to the segment [-h/2,h/2]

  Stress(0,0) = Elascoef*overhalfh*sxxs1;  //  we put directly here the constant E'/(4Pi)
  Stress(0,1) = Elascoef*overhalfh*sxys1;
  Stress(0,2) = Elascoef*overhalfh*syys1;
  Stress(0,3) = Elascoef*overhalfh*syyn1;

  Stress(1,0) = Elascoef*overhalfh*sxxs2;
  Stress(1,1) = Elascoef*overhalfh*sxys2;
  Stress(1,2) = Elascoef*overhalfh*syys2;
  Stress(1,3) = Elascoef*overhalfh*syyn2;

  return Stress;
}


//    Function to get the normal and shear stress at a point on a surface (with given normal and shear vector) induced by
//a linear DD segment of size h  (the linear DD segment is assumed to be the unit element along the x axis [-h/2,h/2]
//Material of Plane-strain Young's modulus Ep
void NormalShearStressKernel_LinearDD(il::StaticArray2D<double,2,4>& St,const il::StaticArray<double,2> xe, const double & h,const il::StaticArray<double,2> s, const il::StaticArray<double,2> n,
                                      const double  Ep ){
  // St:: returned stress
  // xe (x,y):: collocation point to compute the stress in the reference frame of the unit element
  // h :: elt size
  // s :: tangent vector at xe on which to project to obtain the shear stress
  // n:: nornal vector at xe on which to project to obtain the normal stress

  double n1n1,n2n2,n1s1,n2s2,n1n2,n1s2pn2s1 ;
  il::StaticArray2D<double,2,4>  st ;

  n1n1=n[0]*n[0];
  n2n2=n[1]*n[1];
  n1n2=n[0]*n[1];
  n1s1=n[0]*s[0];
  n2s2=n[1]*s[1];

  n1s2pn2s1=n[0]*s[1]+n[1]*s[0];

  st=StressesKernelLinearDD(h,Ep,xe[0],xe[1]);// columns sxxs, sxys, syys, syyn    (knowing that we sxxn and sxyn are respectively equal to sxys and syys )

  double sh1s,sn1s,sh1n,sn1n,sh2s,sn2s,sh2n,sn2n ;
   //shear stress
  // node 1
  sh1s =n1s1*st(0,0)+n1s2pn2s1*st(0,1)+n2s2*st(0,2); //shear dd
  sh1n =n1s1*st(0,1)+n1s2pn2s1*st(0,2)+n2s2*st(0,3) ; // normal dd
  //node 2
  sh2s =n1s1*st(1,0)+ n1s2pn2s1*st(1,1)+n2s2*st(1,2); //shear dd
  sh2n =n1s1*st(1,1)+n1s2pn2s1*st(1,2)+n2s2*st(1,3) ; // normal dd

  //normal stress
  // node 1
  sn1s =n1n1*st(0,0)+2*n1n2*st(0,1)+n2n2*st(0,2); // shear dd
  sn1n = n1n1*st(0,1)+2*n1n2*st(0,2)+n2n2*st(0,3); // normal dd
  //node 2
  sn2s =n1n1*st(1,0)+2*n1n2*st(1,1)+n2n2*st(1,2); // shear dd
  sn2n =n1n1*st(1,1)+2*n1n2*st(1,2)+n2n2*st(1,3); // normal dd

  // output the desired stress components induced either by normal or shear dd of the 2 nodes of the linear DD segment.
  St(0,0)=sh1s;St(0,1)=sh1n;St(0,2)=sh2s;St(0,3)=sh2n;
  St(1,0)=sn1s;St(1,1)=sn1n;St(1,2)=sn2s;St(1,3)=sn2n;
   // St dimensions 2 * 4
}

