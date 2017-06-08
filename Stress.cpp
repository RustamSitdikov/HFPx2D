//
// Created by DONG LIU on 1/30/17.
//
//
//
//Calculate the stress of a given point (xg,yg)(coordinates in the global system)

//this has been ajusted with the modified source code.

#include "Stress.h"

#include <il/math.h>
#include <il/linear_algebra.h>


//(We put these two functions take_submatrix() and set_submatrix() in the headfile
// of the AssemblyDDM.h)

namespace hfp2d {

    il::Array2D<double> stresscalculation(Mesh mesh, il::Array2D<int> id, int p,
                                          const double Ep, double xg,
                                          double yg) {
        // mesh:: the Mesh object
        // id :: the DOF handle
        // p :: the interpolation order
        // Ep :: the Plane Strain Young's modulus
        IL_EXPECT_FAST(id.size(0) == mesh.nelts());
        IL_EXPECT_FAST(id.size(1) == 2 * (p + 1));

        il::Array2D<double> kmat{4, 2 * (p + 1) * id.size(0), 0};
        il::Array<int> dofe{2 * (p + 1), 0};
        il::StaticArray2D<double, 2, 2> R;
        il::Array2D<double> xe{2, 2, 0.};
        //il::StaticArray<double, 2> sec;//no need to be static variable
        SegmentCharacteristic mysege;
        il::StaticArray<double, 2> x_e;//no need to be static variable, but constrained by line 62, il::dot()
        il::StaticArray<double, 2> xe_e;
        il::StaticArray<double, 2> vec01, vec10;
        vec01[0] = 0;
        vec01[1] = 1;
        vec10[0] = 1;
        vec10[1] = 0;
        il::StaticArray2D<double, 2, 4> st1;
        il::StaticArray2D<double, 2, 4> st2;

        for (int e = 0; e < mesh.nelts(); ++e) { // loop on all  elements

            //take_submatrix(xe, mesh.connectivity(e, 0), mesh.connectivity(e, 1), 0, 1,
            //               mesh.coor());
            // take the coordinates of element e from the mesh object
            mysege = get_segment_DD_characteristic(mesh, e, p);
            // get the segment characteristic.
            R = rotation_matrix_2D(mysege.theta);
            //Rotation matrix of the element w.r. to x-axis.

            for (int i = 0;
                 i < 2 * (p + 1); ++i) { // vector of dof id of the element e
                dofe[i] = id(e, i);
            };


            x_e[0] = xg - mysege.Xmid[0];
            x_e[1] = yg - mysege.Xmid[1];

            xe_e = il::dot(R, x_e);

            st1=normal_shear_stress_kernel_dp1_dd(xe_e, mysege.size, vec01,
                                             vec10,
                                             Ep);
            st2=normal_shear_stress_kernel_dp1_dd( xe_e, mysege.size, vec10,
                                             vec01,
                                             Ep);
            for (int ic = 0; ic < 4; ++ic) {
                for (int j = 0; j < 4; ++j) {
                    if (ic < 2) {
                        kmat(ic, dofe[j]) = st1(ic, j);//*Ep/M_PI/4;
                    } else {
                        kmat(ic, dofe[j]) = st2(ic - 2, j);//*Ep/M_PI/4;
                    }
                }
            }
        }
        return kmat;//change this on Februray the 7th, to add the coefficient Ep/pi/4 needs to be verified because in Normalshear funciton there's already the effect of Ep
    }


    il::Array<double> stressoutput(Mesh mesh, il::Array2D<int> id, int p,
                                   const double Ep, double xg, double yg,
                                   il::Array<double> width) {
        il::Array2D<double> kmat;
        kmat = stresscalculation(mesh, id, p, Ep, xg, yg);
        il::Array<double> stressb;
        stressb = il::dot(kmat, width);//(sxy,sxx,syx,syy)
        il::Array<double> stress{3, 0};
        stress[0] = stressb[1];//sxx
        stress[1] = stressb[3];//syy
        stress[2] = stressb[0];//sxy
        return stress;
    }
}