//
// Created by victor on 10/02/2020.
//
#include <iostream>

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzvec.h"

#include "MixedHdiv.h"

//Constants
const int porder = 1;
const int dim = 2;
bool debug = true;

int main(){
    TPZVec<int> nx(2,0), hx(2,0);
    nx[0] = nx[1] = 8; hx[0] = hx[1] = 1;
    TPZVec<REAL> x0 ={0,0,0},x1 = {1,1,0};

    DarcyMixedHdiv MixedSimulator(nx,hx,dim,porder,debug);
    MixedSimulator.SolveMixedHdiv();
}