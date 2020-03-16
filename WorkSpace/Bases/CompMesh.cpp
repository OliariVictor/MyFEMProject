//
// Created by victor on 19/02/2020.
//

#include "CompMesh.h"

Computational::Computational(){
    fdim = 2;
    fmatId = 1;

    fmatBCbottom = -1;
    fmatBCright = -2;
    fmatBCtop = -3;
    fmatBCleft = -4;

    fdirichlet = 0;
    fneumann = 1;
}

TPZCompMesh* Computational::CMesh_m(TPZGeoMesh *gmesh, int pOrder, TLaplaceExample1 &Laplace,TPZFMatrix<REAL>  K, TPZFMatrix<REAL>  invK){

    //Computational Mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh -> SetDefaultOrder(pOrder);
    cmesh -> SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    //Material Creation
    TPZMixedPoisson *material = new TPZMixedPoisson(fmatId,fdim); //Using standard PermealityTensor = Identity.
    if( K(0,0) != 1. || K(1,1) !=1. || K(2,2) != 1.){
        material->SetPermeabilityTensor(K,invK);
    }

    material->SetForcingFunction(Laplace.ForcingFunction());
    material->SetForcingFunctionExact(Laplace.Exact());
    cmesh->InsertMaterialObject(material);


    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);

    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0

    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbottom, fdirichlet , val1, val2); //Cria material que implementa a condição de contorno inferior
    BCond0->SetForcingFunction(Laplace.Exact());
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCright, fneumann , val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunction(Laplace.Exact());
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunction(Laplace.Exact());
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunction(Laplace.Exact());
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    cmesh->AutoBuild();
    cmesh->ExpandSolution();

    return cmesh;
}

TPZCompMesh* Computational::CMeshProperFunc_m(TPZGeoMesh *gmesh, int pOrder,void (*f_source)(const TPZVec<REAL> &x, TPZVec<STATE> &val),void (*Sol_Exact)(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf)){

    //Computational Mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh -> SetDefaultOrder(pOrder);
    cmesh -> SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    //Material Creation
    TPZMixedPoisson *material = new TPZMixedPoisson(fmatId,fdim); //Using standard PermealityTensor = Identity.

    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE>(f_source,5);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Sol_Exact, 5);

    material->SetForcingFunction(fp);
    material->SetForcingFunctionExact(solp);
    cmesh->InsertMaterialObject(material);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);

    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0

    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbottom, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    BCond0->SetForcingFunction(Sol_Exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCright, fneumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunction(Sol_Exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunction(Sol_Exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunction(Sol_Exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond3); //Insere material na malha

    cmesh->AutoBuild();
    cmesh->ExpandSolution();

    return cmesh;
}

TPZCompMesh* Computational::CMesh_flux(TPZGeoMesh *gmesh, int pOrder) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);

    cmesh->SetAllCreateFunctionsHDiv();

    TPZNullMaterial *material = new TPZNullMaterial(fmatId);
    material->SetDimension(fdim);
    cmesh->InsertMaterialObject(material);

    //Create Boundary conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    TPZMaterial *BCond0 = material->CreateBC(material,fmatBCbottom,fdirichlet,val1,val2);
    cmesh->InsertMaterialObject(BCond0);

    TPZMaterial *BCond1 = material->CreateBC(material,fmatBCright,fneumann,val1,val2);
    cmesh->InsertMaterialObject(BCond1);

    TPZMaterial *BCond2 = material->CreateBC(material,fmatBCtop,fdirichlet,val1,val2);
    cmesh->InsertMaterialObject(BCond2);

    TPZMaterial *BCond3 = material->CreateBC(material,fmatBCleft,fdirichlet,val1,val2);
    cmesh->InsertMaterialObject(BCond3);

    cmesh->AutoBuild();
    cmesh->ExpandSolution();

    return cmesh;
}

TPZCompMesh* Computational::CMesh_p(TPZGeoMesh *gmesh, int pOrder){
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(fdim);

    cmesh->SetAllCreateFunctionsContinuous(); //H1 functions
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial *material = new TPZNullMaterial(fmatId); material->SetDimension(fdim);
    cmesh->InsertMaterialObject(material);

    TPZFMatrix<STATE> val1(2,2,0.), val2(2,1,0.);

    cmesh->AutoBuild();
    cmesh->ExpandSolution();

    TPZAdmChunkVector<TPZConnect> &nodeIt = cmesh->ConnectVec();
    for(auto &nodo : nodeIt){
        nodo.SetLagrangeMultiplier(1);
    }
    return cmesh;
}