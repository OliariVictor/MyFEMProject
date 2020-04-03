//
// Created by victor on 10/02/2020.
//
#include <iostream>

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"
#include "pzfmatrix.h"

#include "MixedHdiv.h"
#include "Loads.h"
#include "TPZAnalyticSolution.h"

//Constants
const int porder = 1;
const int dim = 2;
bool debug = true;

int main(){

    TPZFNMatrix<9,REAL> K = {{2,0,0},{0,5,0},{0,0,1}}, invK = {{0.5,0,0},{0,0.2,0},{0,0,1}};

    //TLaplaceExample1 Lp;
    TLaplaceExample1 Lp(K,invK);
    Lp.fExact = TLaplaceExample1::ESinSin;
    Lp.fSignConvention= 1;
    bool postproc = false;  // PostProc computational cost is too high

    //If enabled compute an Hdiv approximation for the specified mesh
    if(0) {
        TPZVec<int> nx(2,0), hx(2,0);
        nx[0] = nx[1] = 2; hx[0] = hx[1] = 1;

        DarcyMixedHdiv MixedSimulator(nx, hx, dim, porder, debug);
        MixedSimulator.SolveMixedHdiv(Lp);
    }

    if(1){
        int numRefinements = 3;

        DarcyMixedHdiv errorRate(porder);
        //errorRate.ErrorRate(numRefinements,Lp);
        errorRate.ErrorRate(numRefinements,Lp,K,invK,postproc);
    }
}