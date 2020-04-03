//
// Created by victor on 19/02/2020.
//

#ifndef HDIVEX_DARCYHDIV_H
#define HDIVEX_DARCYHDIV_H

#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzvisualmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZAnalyticSolution.h"

#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZHybridizeHDiv.h"

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
    DarcyMixedHdiv(int pOrder);

    DarcyMixedHdiv(TPZVec<int> nx,TPZVec<int> hx,int dim,int pOrder,bool debug);

    ~DarcyMixedHdiv();

    //Solve a Mixed Hdiv approximation for a mesh previously defined
    void SolveMixedHdiv(TLaplaceExample1 &Laplace);

    void SolveMixedHdivProperFunc(void (*f_source)(const TPZVec<REAL> &x, TPZVec<STATE> &val),void (*Sol_Exact)(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf));

    void SolveMixedProblem(TPZCompMesh *mesh, bool optBW, TLaplaceExample1 &Laplace, TPZVec<REAL> &errorsLog, TPZVec<REAL> &rt, TPZGeoMesh* gmesh, REAL &hLog, REAL &h,bool &last);

    void SolveHybridProblem(TPZCompMesh *mesh, bool optBW, TLaplaceExample1 &Laplace, TPZVec<REAL> &errorsLog, TPZVec<REAL> &rt, TPZGeoMesh* gmesh, REAL &hLog, REAL &h, bool &last);

    void PostProcess(string fileName, int dim, TPZAnalysis &analysis);

    void ErrorRate(int numRefinement,TLaplaceExample1 &,TPZFMatrix<REAL>  K ={{1,0,0},{0,1,0},{0,0,1}}, TPZFMatrix<REAL>  invK ={{1,0,0},{0,1,0},{0,0,1}}, bool postproc = false);

    void ErrorRateProperFunc(int numRefinement,void (*f_source)(const TPZVec<REAL> &x, TPZVec<STATE> &val),void (*Sol_Exact)(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf));
};

#endif //HDIVEX_DARCYHDIV_H
