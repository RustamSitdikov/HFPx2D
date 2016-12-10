//
// Created by Brice Lecampion on 10.12.16.
//

#include "Mesh.h"


void Mesh::set_values(il::Array2D<double> xy, il::Array2D<int> ien)
{
  IL_ASSERT(xy.size(1) ==2); //check array dimensions ?
  IL_ASSERT(ien.size(1) ==2);//check array dimensions ??? -> this is only for 1D mesh so far

  Coor = xy;   // list of coordinates of points in the mesh
  conn = ien;  //  connectivity array -
}
// needs to add function to add one or more elements ... (needs to have active and passive elements to track active/passive fractures etc.)
//
/// // could provide a default constructor for a straight fracture ?

