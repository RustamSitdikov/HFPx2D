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


namespace hfp2d{

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

//

il::Array2D<il::int_t> GetEdgesSharing2(hfp2d::Mesh &mesh ){
  // get element sharing common edges

  // loop over all nodes of the mesh
  // find corresponding elements
  // if n of elt sharing that node ==2 -> put in array...
  //

  // maximum number of vertex with 2 elements is
  // this will not work for the case of more than 2 elements sharing the node.
  // we don t care of that case for now.


  // format is col1: element1, col2: element2   (note we don t store the corresponding node here....)

  il::Array2D<il::int_t> edge{mesh.numberOfNodes(),2,0.};

  il::StaticArray<il::int_t,2> temp;

  int j=0; int k=0;

  for (il::int_t i=0;i<mesh.numberOfNodes();i++){

    j=0;
//    temp[j]=i;j++;
     for (il::int_t e=0;e<mesh.numberOfElements();e++){

        if (mesh.connectivity(e,0)==i) {
          temp[j]=e;
          j++;
        }
       else if (mesh.connectivity(e,1)==i){
          temp[j]=e;
          j++;
        };

       if (j==2) {
         edge(k,0) = temp[0];
         edge(k, 1) = temp[1];
//         edge(k, 2) = temp[2];
//         temp[0]=0;temp(0,1)=0;
         j = 0;
         k++;
         break;
       };

     }

    j=0;

  }

  il::Array2D<il::int_t> outputedge{k,2,0.};

  for (il::int_t i=0;i<outputedge.size(0);i++) {

    outputedge(i,0)=edge(i,0);
    outputedge(i,1)=edge(i,1);
//    outputedge(i,2)=edge(i,2);

  }

  return outputedge;
}






// FV Matrix L

// should return a sparse matrix....
// for now , for quick debug assemble a dense ;(

il::Array2D<double> BuildFD_P0(hfp2d::Mesh &mesh, hfp2d::Fluid &fluid, il::Array<double> &hydraulicwidth){
//
//
//
// checks here....


 // loop on edges.....
  il::Array2D<il::int_t> edgecommon=hfp2d::GetEdgesSharing2(mesh);

  il::Array2D<double> L{mesh.numberOfElements(),mesh.numberOfElements(),0.}; // should be a sparse matrix.....

  // compute conductivities.....
  il::Array<double> CurrentCond=EdgeConductivitiesP0Newtonian(edgecommon,hydraulicwidth,fluid);

  il::int_t er,el,ml;
  double hi;
  ml=edgecommon.size(0);

  for (il::int_t i=0; i<edgecommon.size(0);i++ ){

    er=edgecommon(i,0);
    el=edgecommon(i,1);

    hi = (mesh.eltsize(er)+mesh.eltsize(el))/2.;

    L(er,er)+=CurrentCond[i]/hi;
    L(el,el)+=CurrentCond[i]/hi;
    L(er,el)=-CurrentCond[i]/hi;
    L(el,er)=-CurrentCond[i]/hi;

  };

  return L;

 }


// function for current cell volume -> returning a vector

// function for elt size  -> returning a vector



}