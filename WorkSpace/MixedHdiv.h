//
// Created by victor on 19/02/2020.
//

#ifndef HDIVEX_DARCYHDIV_H
#define HDIVEX_DARCYHDIV_H

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzvisualmatrix.h"

#include "TPZVTKGeoMesh.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"

#include "GeoMesh.h"
#include "CompMesh.h"
#include "Loads.h"

class DarcyMixedHdiv{
private:
    TPZVec<int> fnx;
    TPZVec<int> fhx;
    int fPOrder;
    int fdim;
    bool fdebug;
public:
    DarcyMixedHdiv(TPZVec<int> nx,TPZVec<int> hx,int dim,int pOrder,bool debug);

    ~DarcyMixedHdiv();

    void SolveMixedHdiv();
};

#endif //HDIVEX_DARCYHDIV_H
