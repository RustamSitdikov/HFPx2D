//
// This file is part of HFPx2DUnitTest.
//
// Created by Brice Lecampion on 06.10.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#include "ReynoldsP0.h"

#include "il/Array.h"

#include "src/core/Mesh.h"
#include "src/core/Fluid.h"


namespace hfp2D{

il::Array<double>  EdgeConductivitiesP0Newtonian(il::Array2D<il::int_t> &edgeAdj,il::Array<double> &hydraulic_width,hfp2d::Fluid &fluid){
  // input should be an array of int : row edges, column the 2 adjacent elements . (and we assume that we element# = width(dof#) for simplicity
  //  hydraulic_width array with hydraulic width in each elements.

  double poiseuille_ct= 1./(12.*fluid.fluidViscosity()) ;
  double cr,cl;

  il::Array<double> Cond{edgeAdj.size(0),0};
  // loop on the inner edges

  // harmonic mean   2*w_l^3*w_r^3/(w_l^3 + w_r^3)
  for (il::int_t i=0;i<edgeAdj.size(0);i++ ){
    cl=pow(hydraulic_width[edgeAdj(i,0)],3.);
    cr=pow(hydraulic_width[edgeAdj(i,1)],3.);

    Cond[i]=poiseuille_ct*2.*(cr*cl)/(cr+cl);

  }

  return Cond;

};


// FV Matrix L

// should return a sparse matrix....
// for now , for quick debug assemble a dense ;(

il::Array2D<double> BuildFD_P0(hfp2d::Mesh &mesh, il::Array2D<il::int_t> &edgeAdj, il::Array<double> &Cond ){



 }






}