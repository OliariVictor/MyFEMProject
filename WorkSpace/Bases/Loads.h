//
// Created by victor on 19/02/2020.
//

#ifndef HDIVEX_LOADS_H
#define HDIVEX_LOADS_H

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "tpzautopointer.h"

namespace Loads {

    void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f);

    void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
}
#endif //HDIVEX_LOADS_H
