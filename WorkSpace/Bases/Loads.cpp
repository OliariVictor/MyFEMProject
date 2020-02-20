//
// Created by victor on 19/02/2020.
//

#include "Loads.h"


void Loads::Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){

    dsol.Resize(3,1);
    sol.Resize(3);
    const REAL Pi=M_PI;

    REAL xv = x[0];
    REAL yv = x[1];

    STATE sigma_x = -2.*Pi*cos(2.*Pi*xv)*sin(2.*Pi*yv);
    STATE sigma_y = -2.*Pi*cos(2.*Pi*yv)*sin(2.*Pi*xv);
    STATE pressure= sin(2.*Pi*xv)*sin(2.*Pi*yv);

    sol[0]=pressure;

    // Gradiente pressão

    dsol(0,0)= -sigma_x;
    dsol(1,0)= -sigma_y;
}

void Loads::F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f){

    f.resize(1);
    const REAL Pi=M_PI;

    REAL xv = x[0];
    REAL yv = x[1];

    STATE f_x = -8.0*Pi*Pi*sin(2.0*Pi*xv)*sin(2.0*Pi*yv);

    f[0] = f_x;
}