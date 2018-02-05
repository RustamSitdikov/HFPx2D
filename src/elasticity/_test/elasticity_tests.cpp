//
// This file is part of HFPx2D.
//
// Created by Federico Ciardo on 05.02.18.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2018.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include <iostream>

#include <gtest/gtest.h>

#include <il/Array2D.h>
#include <il/linear_algebra.h>
#include <il/norm.h>

#include <src/elasticity/AssemblyDDM.h>
#include <src/elasticity/PlaneStrainInfinite.h>
#include <src/elasticity/Simplified3D.h>
#include <src/core/Mesh.h>

//--------------------------------------------------------------------------
// P1 Plane-Strain

/// TEST 1

TEST(P1, two_segs_45_a1) {
    // two segments at 90 degree from one another
    // oriented 45 degree from the axis e_x
    //
    //        /\        //

    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 1.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K = hfp2d::basic_assembly(
            mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

    // we compare the results of the assembly w.t the mma code
    il::Array2D<double> Kmma{8, 8, 0.};
    ;
    Kmma(0, 0) = 0.483423;
    Kmma(0, 1) = -1.20657e-16;
    Kmma(0, 2) = -0.0332652;
    Kmma(0, 3) = 2.07014e-17;
    Kmma(0, 4) = 0.00824514;
    Kmma(0, 5) = -0.0381383;
    Kmma(0, 6) = -0.0133572;
    Kmma(0, 7) = -0.0321492;

    Kmma(1, 0) = -2.27998e-16;
    Kmma(1, 1) = 0.483423;
    Kmma(1, 2) = 2.80878e-17;
    Kmma(1, 3) = -0.0332652;
    Kmma(1, 4) = -0.00355156;
    Kmma(1, 5) = -0.00824514;
    Kmma(1, 6) = 0.00954068;
    Kmma(1, 7) = 0.0133572;

    Kmma(2, 0) = -0.0332652;
    Kmma(2, 1) = 1.4713e-33;
    Kmma(2, 2) = 0.483423;
    Kmma(2, 3) = 0.;
    Kmma(2, 4) = -0.0536082;
    Kmma(2, 5) = -0.376167;
    Kmma(2, 6) = 0.000833234;
    Kmma(2, 7) = -0.0157962;

    Kmma(3, 0) = 7.38637e-18;
    Kmma(3, 1) = -0.0332652;
    Kmma(3, 2) = -1.07342e-16;
    Kmma(3, 3) = 0.483423;
    Kmma(3, 4) = 0.23189;
    Kmma(3, 5) = 0.0536082;
    Kmma(3, 6) = 0.128481;
    Kmma(3, 7) = -0.000833234;

    Kmma(4, 0) = 0.000833234;
    Kmma(4, 1) = 0.0157962;
    Kmma(4, 2) = -0.0536082;
    Kmma(4, 3) = 0.376167;
    Kmma(4, 4) = 0.483423;
    Kmma(4, 5) = 0.;
    Kmma(4, 6) = -0.0332652;
    Kmma(4, 7) = 0.;

    Kmma(5, 0) = -0.128481;
    Kmma(5, 1) = -0.000833234;
    Kmma(5, 2) = -0.23189;
    Kmma(5, 3) = 0.0536082;
    Kmma(5, 4) = 1.07342e-16;
    Kmma(5, 5) = 0.483423;
    Kmma(5, 6) = -7.38637e-18;
    Kmma(5, 7) = -0.0332652;

    Kmma(6, 0) = -0.0133572;
    Kmma(6, 1) = 0.0321492;
    Kmma(6, 2) = 0.00824514;
    Kmma(6, 3) = 0.0381383;
    Kmma(6, 4) = -0.0332652;
    Kmma(6, 5) = -2.07014e-17;
    Kmma(6, 6) = 0.483423;
    Kmma(6, 7) = 1.20657e-16;

    Kmma(7, 0) = -0.00954068;
    Kmma(7, 1) = 0.0133572;
    Kmma(7, 2) = 0.00355156;
    Kmma(7, 3) = -0.00824514;
    Kmma(7, 4) = -2.80878e-17;
    Kmma(7, 5) = -0.0332652;
    Kmma(7, 6) = 2.27998e-16;
    Kmma(7, 7) = 0.483423;

    std::cout << K(0, 0) << "\n";
    std::cout << K(0, 2) << "\n";

    double my_sum = 0.;

    for (il::int_t j = 0.; j < K.size(1); j++) {
        for (il::int_t i = 0; i < K.size(0); i++) {
            my_sum += abs(K(i, j) - Kmma(i, j));
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_NEAR(my_sum, 0., 1.e-5);
}

/// TEST 2

TEST(P1, two_segs_90_45_a1) {
    // two segments at 90 degree from one another
    // oriented 45 degree from the axis e_x
    //
    //        /\        //

    // compare with the case of 2 segmenets at 90 from one another
    // but oriented along e_x and e_y
    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 1.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh45(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K45 = hfp2d::basic_assembly(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

    il::Array2D<double> xy90{3, 2, 0.};
    xy90(0, 0) = 0.;
    xy90(0, 1) = 0.;
    xy90(1, 0) = sqrt(2.);
    xy90(1, 1) = 0.;
    xy90(2, 0) = sqrt(2.);
    xy90(2, 1) = -sqrt(2.);

    hfp2d::Mesh mesh90(xy90, ien, 1);

    il::Array2D<double> K90 = hfp2d::basic_assembly(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

    double my_sum = 0.;

    for (il::int_t j = 0.; j < K45.size(1); j++) {
        for (il::int_t i = 0; i < K45.size(0); i++) {
            my_sum += abs(K45(i, j) - K90(i, j));
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_TRUE(my_sum == 0.);
}

/// TEST 3

TEST(P1, two_adjacent_segs) {
    // two adjacents straight segments
    // just one DD is mobilised
    //        ._._.        //

    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 0.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K = hfp2d::basic_assembly(
            mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd, 0.);

    // we compare the results of the assembly w.t the mma code
    il::Array2D<double> Kmma{8, 8, 0.};
    ;
    Kmma(0, 0) = -0.683664;
    Kmma(0, 1) = 0.;
    Kmma(0, 2) = 0.0470442;
    Kmma(0, 3) = 0.;
    Kmma(0, 4) = 0.0943392;
    Kmma(0, 5) = 0.;
    Kmma(0, 6) = 0.0187761;
    Kmma(0, 7) = 0.;

    Kmma(1, 0) = 0.;
    Kmma(1, 1) = -0.683664;
    Kmma(1, 2) = 0.;
    Kmma(1, 3) = 0.0470442;
    Kmma(1, 4) = 0.;
    Kmma(1, 5) = 0.0943392;
    Kmma(1, 6) = 0.;
    Kmma(1, 7) = 0.0187761;

    Kmma(2, 0) = 0.0470442;
    Kmma(2, 1) = 0.;
    Kmma(2, 2) = -0.683664;
    Kmma(2, 3) = 0.;
    Kmma(2, 4) = 0.379637;
    Kmma(2, 5) = 0.;
    Kmma(2, 6) = 0.0315223;
    Kmma(2, 7) = 0.;

    Kmma(3, 0) = 0.;
    Kmma(3, 1) = 0.0470442;
    Kmma(3, 2) = 0.;
    Kmma(3, 3) = -0.683664;
    Kmma(3, 4) = 0.;
    Kmma(3, 5) = 0.379637;
    Kmma(3, 6) = 0.;
    Kmma(3, 7) = 0.0315223;

    Kmma(4, 0) = 0.0315223;
    Kmma(4, 1) = 0.;
    Kmma(4, 2) = 0.379637;
    Kmma(4, 3) = 0.;
    Kmma(4, 4) = -0.683664;
    Kmma(4, 5) = 0.;
    Kmma(4, 6) = 0.0470442;
    Kmma(4, 7) = 0.;

    Kmma(5, 0) = 0.;
    Kmma(5, 1) = 0.0315223;
    Kmma(5, 2) = 0.;
    Kmma(5, 3) = 0.379637;
    Kmma(5, 4) = 0.;
    Kmma(5, 5) = -0.683664;
    Kmma(5, 6) = 0.;
    Kmma(5, 7) = 0.0470442;

    Kmma(6, 0) = 0.0187761;
    Kmma(6, 1) = 0.;
    Kmma(6, 2) = 0.0943392;
    Kmma(6, 3) = 0.;
    Kmma(6, 4) = 0.0470442;
    Kmma(6, 5) = 0.;
    Kmma(6, 6) = -0.683664;
    Kmma(6, 7) = 0.;

    Kmma(7, 0) = 0.;
    Kmma(7, 1) = 0.0187761;
    Kmma(7, 2) = 0.;
    Kmma(7, 3) = 0.0943392;
    Kmma(7, 4) = 0.;
    Kmma(7, 5) = 0.0470442;
    Kmma(7, 6) = 0.;
    Kmma(7, 7) = -0.683664;

    double my_sum = 0.;

    for (il::int_t j = 0; j < K.size(1); j++) {
        for (il::int_t i = 0; i < K.size(0); i++) {
            my_sum += K(i, j) + Kmma(i, j);
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_NEAR(my_sum, 0., 1.e-5);
}

//--------------------------------------------------------------------------
// P0 simplified 3D

/// TEST 4

TEST(P0, two_segs_90_45_a1) {
    // two segments at 90 degree from one another
    // oriented 45 degree from the axis e_x
    //
    //        /\        //

    // compare with the case of 2 segmenets at 90 from one another
    // but oriented along e_x and e_y
    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 1.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh45(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K45 = hfp2d::basic_assembly(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, 10.);

    il::Array2D<double> xy90{3, 2, 0.};
    xy90(0, 0) = 0.;
    xy90(0, 1) = 0.;
    xy90(1, 0) = sqrt(2.);
    xy90(1, 1) = 0.;
    xy90(2, 0) = sqrt(2.);
    xy90(2, 1) = -sqrt(2.);

    hfp2d::Mesh mesh90(xy90, ien, 1);

    il::Array2D<double> K90 = hfp2d::basic_assembly(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, 10.);

    double my_sum = 0.;

    for (il::int_t j = 0.; j < K45.size(1); j++) {
        for (il::int_t i = 0; i < K45.size(0); i++) {
            my_sum += abs(K45(i, j) - K90(i, j));
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_TRUE(my_sum == 0.);
}

/// TEST 5

TEST(P1_nodal, two_segs_45_a1) {
    // two segments at 90 degree from one another
    // oriented 45 degree from the axis e_x
    //
    //        /\        //

    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 1.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K = hfp2d::basic_assembly_nodal(
            mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd_nodal, 0.);

    // we compare the results of the assembly w.t the mma code
    il::Array2D<double> Kmma{8, 8, 0.};
    ;
    Kmma(0, 0) = 0.483423;
    Kmma(0, 1) = -1.20657e-16;
    Kmma(0, 2) = -0.0332652;
    Kmma(0, 3) = 2.07014e-17;
    Kmma(0, 4) = 0.00824514;
    Kmma(0, 5) = -0.0381383;
    Kmma(0, 6) = -0.0133572;
    Kmma(0, 7) = -0.0321492;

    Kmma(1, 0) = -2.27998e-16;
    Kmma(1, 1) = 0.483423;
    Kmma(1, 2) = 2.80878e-17;
    Kmma(1, 3) = -0.0332652;
    Kmma(1, 4) = -0.00355156;
    Kmma(1, 5) = -0.00824514;
    Kmma(1, 6) = 0.00954068;
    Kmma(1, 7) = 0.0133572;

    Kmma(2, 0) = -0.0332652;
    Kmma(2, 1) = 1.4713e-33;
    Kmma(2, 2) = 0.483423;
    Kmma(2, 3) = 0.;
    Kmma(2, 4) = -0.0536082;
    Kmma(2, 5) = -0.376167;
    Kmma(2, 6) = 0.000833234;
    Kmma(2, 7) = -0.0157962;

    Kmma(3, 0) = 7.38637e-18;
    Kmma(3, 1) = -0.0332652;
    Kmma(3, 2) = -1.07342e-16;
    Kmma(3, 3) = 0.483423;
    Kmma(3, 4) = 0.23189;
    Kmma(3, 5) = 0.0536082;
    Kmma(3, 6) = 0.128481;
    Kmma(3, 7) = -0.000833234;

    Kmma(4, 0) = 0.000833234;
    Kmma(4, 1) = 0.0157962;
    Kmma(4, 2) = -0.0536082;
    Kmma(4, 3) = 0.376167;
    Kmma(4, 4) = 0.483423;
    Kmma(4, 5) = 0.;
    Kmma(4, 6) = -0.0332652;
    Kmma(4, 7) = 0.;

    Kmma(5, 0) = -0.128481;
    Kmma(5, 1) = -0.000833234;
    Kmma(5, 2) = -0.23189;
    Kmma(5, 3) = 0.0536082;
    Kmma(5, 4) = 1.07342e-16;
    Kmma(5, 5) = 0.483423;
    Kmma(5, 6) = -7.38637e-18;
    Kmma(5, 7) = -0.0332652;

    Kmma(6, 0) = -0.0133572;
    Kmma(6, 1) = 0.0321492;
    Kmma(6, 2) = 0.00824514;
    Kmma(6, 3) = 0.0381383;
    Kmma(6, 4) = -0.0332652;
    Kmma(6, 5) = -2.07014e-17;
    Kmma(6, 6) = 0.483423;
    Kmma(6, 7) = 1.20657e-16;

    Kmma(7, 0) = -0.00954068;
    Kmma(7, 1) = 0.0133572;
    Kmma(7, 2) = 0.00355156;
    Kmma(7, 3) = -0.00824514;
    Kmma(7, 4) = -2.80878e-17;
    Kmma(7, 5) = -0.0332652;
    Kmma(7, 6) = 2.27998e-16;
    Kmma(7, 7) = 0.483423;

    std::cout << K(0, 0) << "\n";
    std::cout << K(0, 2) << "\n";

    double my_sum = 0.;

    for (il::int_t j = 0.; j < K.size(1); j++) {
        for (il::int_t i = 0; i < K.size(0); i++) {
            my_sum += abs(K(i, j) - Kmma(i, j));
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_NEAR(my_sum, 0., 1.e-5);
}

/// TEST 6

TEST(P1_nodal, two_segs_90_45_a1) {
    // two segments at 90 degree from one another
    // oriented 45 degree from the axis e_x
    //
    //        /\        //

    // compare with the case of 2 segmenets at 90 from one another
    // but oriented along e_x and e_y
    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 1.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh45(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K45 = hfp2d::basic_assembly_nodal(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd_nodal, 0.);

    il::Array2D<double> xy90{3, 2, 0.};
    xy90(0, 0) = 0.;
    xy90(0, 1) = 0.;
    xy90(1, 0) = sqrt(2.);
    xy90(1, 1) = 0.;
    xy90(2, 0) = sqrt(2.);
    xy90(2, 1) = -sqrt(2.);

    hfp2d::Mesh mesh90(xy90, ien, 1);

    il::Array2D<double> K90 = hfp2d::basic_assembly_nodal(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd_nodal, 0.);

    double my_sum = 0.;

    for (il::int_t j = 0.; j < K45.size(1); j++) {
        for (il::int_t i = 0; i < K45.size(0); i++) {
            my_sum += abs(K45(i, j) - K90(i, j));
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_TRUE(my_sum == 0.);
}

/// TEST 7

TEST(P1_nodal, two_adjacent_segs) {
    // two adjacents straight segments
    // just one DD is mobilised
    //        ._._.        //

    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 0.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K = hfp2d::basic_assembly_nodal(
            mesh, myelas, hfp2d::normal_shear_stress_kernel_dp1_dd_nodal, 0.);

    // we compare the results of the assembly w.t the mma code
    il::Array2D<double> Kmma{8, 8, 0.};
    ;
    Kmma(0, 0) = -0.683664;
    Kmma(0, 1) = 0.;
    Kmma(0, 2) = 0.0470442;
    Kmma(0, 3) = 0.;
    Kmma(0, 4) = 0.0943392;
    Kmma(0, 5) = 0.;
    Kmma(0, 6) = 0.0187761;
    Kmma(0, 7) = 0.;

    Kmma(1, 0) = 0.;
    Kmma(1, 1) = -0.683664;
    Kmma(1, 2) = 0.;
    Kmma(1, 3) = 0.0470442;
    Kmma(1, 4) = 0.;
    Kmma(1, 5) = 0.0943392;
    Kmma(1, 6) = 0.;
    Kmma(1, 7) = 0.0187761;

    Kmma(2, 0) = 0.0470442;
    Kmma(2, 1) = 0.;
    Kmma(2, 2) = -0.683664;
    Kmma(2, 3) = 0.;
    Kmma(2, 4) = 0.379637;
    Kmma(2, 5) = 0.;
    Kmma(2, 6) = 0.0315223;
    Kmma(2, 7) = 0.;

    Kmma(3, 0) = 0.;
    Kmma(3, 1) = 0.0470442;
    Kmma(3, 2) = 0.;
    Kmma(3, 3) = -0.683664;
    Kmma(3, 4) = 0.;
    Kmma(3, 5) = 0.379637;
    Kmma(3, 6) = 0.;
    Kmma(3, 7) = 0.0315223;

    Kmma(4, 0) = 0.0315223;
    Kmma(4, 1) = 0.;
    Kmma(4, 2) = 0.379637;
    Kmma(4, 3) = 0.;
    Kmma(4, 4) = -0.683664;
    Kmma(4, 5) = 0.;
    Kmma(4, 6) = 0.0470442;
    Kmma(4, 7) = 0.;

    Kmma(5, 0) = 0.;
    Kmma(5, 1) = 0.0315223;
    Kmma(5, 2) = 0.;
    Kmma(5, 3) = 0.379637;
    Kmma(5, 4) = 0.;
    Kmma(5, 5) = -0.683664;
    Kmma(5, 6) = 0.;
    Kmma(5, 7) = 0.0470442;

    Kmma(6, 0) = 0.0187761;
    Kmma(6, 1) = 0.;
    Kmma(6, 2) = 0.0943392;
    Kmma(6, 3) = 0.;
    Kmma(6, 4) = 0.0470442;
    Kmma(6, 5) = 0.;
    Kmma(6, 6) = -0.683664;
    Kmma(6, 7) = 0.;

    Kmma(7, 0) = 0.;
    Kmma(7, 1) = 0.0187761;
    Kmma(7, 2) = 0.;
    Kmma(7, 3) = 0.0943392;
    Kmma(7, 4) = 0.;
    Kmma(7, 5) = 0.0470442;
    Kmma(7, 6) = 0.;
    Kmma(7, 7) = -0.683664;

    double my_sum = 0.;

    for (il::int_t j = 0; j < K.size(1); j++) {
        for (il::int_t i = 0; i < K.size(0); i++) {
            my_sum += K(i, j) + Kmma(i, j);
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_NEAR(my_sum, 0., 1.e-5);
}

//--------------------------------------------------------------------------
// P0 simplified 3D

/// TEST 8

TEST(P0_nodal, two_segs_90_45_a1) {
    // two segments at 90 degree from one another
    // oriented 45 degree from the axis e_x
    //
    //        /\        //

    // compare with the case of 2 segmenets at 90 from one another
    // but oriented along e_x and e_y
    il::Array2D<double> xy{3, 2, 0.};
    xy(0, 0) = -1.;
    xy(0, 1) = 0.;
    xy(1, 0) = 0.;
    xy(1, 1) = 1.;
    xy(2, 0) = 1.;
    xy(2, 1) = 0.;

    il::Array2D<il::int_t> ien{2, 2, 0};
    ien(0, 0) = 0;
    ien(0, 1) = 1;
    ien(1, 0) = 1;
    ien(1, 1) = 2;

    hfp2d::Mesh mesh45(xy, ien, 1);

    hfp2d::ElasticProperties myelas(1., 0.);

    il::Array2D<double> K45 = hfp2d::basic_assembly(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, 10.);

    il::Array2D<double> xy90{3, 2, 0.};
    xy90(0, 0) = 0.;
    xy90(0, 1) = 0.;
    xy90(1, 0) = sqrt(2.);
    xy90(1, 1) = 0.;
    xy90(2, 0) = sqrt(2.);
    xy90(2, 1) = -sqrt(2.);

    hfp2d::Mesh mesh90(xy90, ien, 1);

    il::Array2D<double> K90 = hfp2d::basic_assembly(
            mesh45, myelas, hfp2d::normal_shear_stress_kernel_s3d_dp0_dd, 10.);

    double my_sum = 0.;

    for (il::int_t j = 0.; j < K45.size(1); j++) {
        for (il::int_t i = 0; i < K45.size(0); i++) {
            my_sum += abs(K45(i, j) - K90(i, j));
        }
    }

    std::cout << my_sum << "\n";

    ASSERT_TRUE(my_sum == 0.);
}