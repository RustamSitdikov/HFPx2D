//
// Created by Brice Lecampion on 03.01.17.
//

#ifndef HFPX2D_ASSEMBLYDDM_H
#define HFPX2D_ASSEMBLYDDM_H

#include <il/Array2D.h>
#include "Mesh.h"

void BasicAssembly(il::Array2D<double> &Kmat, Mesh mesh, il::Array2D<int> id, int p , double Ep );


#endif //HFPX2D_ASSEMBLYDDM_H
