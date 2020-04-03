//
// Created by victor on 19/02/2020.
//

#ifndef HDIVEX_COMPMESH_H
#define HDIVEX_COMPMESH_H

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "tpzautopointer.h"
#include "TPZAnalyticSolution.h"

#include "TPZMultiphysicsCompMesh.h"
#include "pzcmesh.h"
#include "mixedpoisson.h"
#include "TPZNullMaterial.h"
#include "pzbndcond.h"

class Computational{
private:

    //General properties
    int fdim;
    int fmatId;
    //Boundary conditions
    int fmatBCbottom;
    int fmatBCright;
    int fmatBCtop;
    int fmatBCleft;
    //Boundary Indices
    int fdirichlet;
    int fneumann;
    //Integration
    int fIntOrder;

public:

    //Constructor
    Computational();

    /**
     * @brief Creates the multiphisics computational mesh of the problem
     * @param *gmesh Pointer to the geometric mesh
     * @param pOrder Integration order of the problem
     */

    TPZMultiphysicsCompMesh *CMesh_m(TPZGeoMesh *gmesh,int pOrder, TLaplaceExample1 &Laplace,TPZFMatrix<REAL>  K ={{1,0,0},{0,1,0},{0,0,1}}, TPZFMatrix<REAL>  invK ={{1,0,0},{0,1,0},{0,0,1}});

    TPZCompMesh* CMeshProperFunc_m(TPZGeoMesh *gmesh, int pOrder,void (*f_source)(const TPZVec<REAL> &x, TPZVec<STATE> &val),void (*Sol_Exact)(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf));

    /**
     * @brief Creates a computational mesh according to the space HDiv
     * @param *gmesh Pointer to the geometric mesh
     * @param pOrder Integration order of the problem
     */
    TPZCompMesh *CMesh_flux(TPZGeoMesh *gmesh, int pOrder);

    /**
     * @brief Creates a computational mesh belonging to a partitioned H1 space
     * @param *gmesh Pointer to geometric mesh
     * @param pOder Integration order of the problem
     */
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int pOrder);
};

#endif //HDIVEX_COMPMESH_H