//
// This file is part of HFPx2D.
//
// Created by Brice Lecampion on 30.08.17.
// Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
// Geo-Energy Laboratory, 2016-2017.  All rights reserved.
// See the LICENSE.TXT file for more details. 
//

#ifndef HFPX2D_SIMPLEELASTIC_H
#define HFPX2D_SIMPLEELASTIC_H


namespace hfp2d {

double SimpleGriffithExampleLinearElement(int nelts);

double SimpleGriffithExampleS3D_P0(int nelts);

double SimpleGriffithExampleLinearElement_AddMesh(int nelts);

double SimpleGriffithExampleS3D_P0_AddMesh(int nelts);

double SimpleGriffithExampleLinearElement_byNodes(int nelts);

}

#endif //HFPX2D_SIMPLEELASTIC_H
