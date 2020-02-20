//
// Created by victor on 18/02/2020.
//

#ifndef HDIVEX_GEOMESH_H
#define HDIVEX_GEOMESH_H

#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzgengrid.h"
#include "TPZVTKGeoMesh.h"
#include "pzcheckgeom.h"

class Geometry{
private:
    //General properties
    int fdim;
    int fmatId;
    //Boundary conditions
    int fmatBCbottom;
    int fmatBCright;
    int fmatBCtop;
    int fmatBCleft;

public:

    //Constructor
    Geometry();

    Geometry(int dim);

/**
 * @brief Hardcoded triangular mesh generator
 * @note This method was created to check the result generated by gengrid
 * @param nx Number of nodes in x direction
 * @param ny Number of nodes in y direction
 * @param hx Width of the mesh
 * @param hy Height of the mesh
 * @return
 */
TPZGeoMesh* gmeshHardWay(int nx, int ny, double hx, double hy);

/**
 * @brief Triangular mesh generated by gengrid
 * @note  The domain was previously partitioned into quadrilaterals and then each quadrilateral was partitioned into two triangles
 * @param nx Number of quadrilaterals in x
 * @param ny Number of quadrilaterals in y direction
 * @param hx Width of the mesh
 * @param hy Height of the mesh
 * @return
 */
TPZGeoMesh* gmeshGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1,int matId);

/**
 * @brief Generates VTK extensions of several grids created by gmeshHardWay and gmeshGenGrid
 * @return
 */
void TestHWvsGG();

};

#endif //HDIVEX_GEOMESH_H