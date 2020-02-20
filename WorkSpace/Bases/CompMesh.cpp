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

TPZCompMesh* Computational::CMesh_m(TPZGeoMesh *gmesh, int pOrder){

    //Computational Mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh -> SetDefaultOrder(pOrder);
    cmesh -> SetDimModel(fdim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();

    //Material Creation
    TPZMixedPoisson *material = new TPZMixedPoisson(fmatId,fdim); //Using standard PermealityTensor = Identity.
    //??? Why one can assign a TPZDummyFunction<> to a TPZAutoPointer<TPZFunction<>>???
    TPZAutoPointer<TPZFunction<STATE> > fp = new TPZDummyFunction<STATE>(Loads::F_source,5);
    TPZAutoPointer<TPZFunction<STATE> > solp = new TPZDummyFunction<STATE> (Loads::Sol_exact, 5);

    material->SetForcingFunction(fp);
    material->SetForcingFunctionExact(solp);
    cmesh->InsertMaterialObject(material);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(2,2,0.), val2(3,1,0.);

    val2(0,0) = 0.0; // vx -> 0
    val2(1,0) = 0.0; // vy -> 0

    TPZMaterial * BCond0 = material->CreateBC(material, fmatBCbottom, fdirichlet, val1, val2); //Cria material que implementa a condição de contorno inferior
    BCond0->SetForcingFunction(Loads::Sol_exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond0); //Insere material na malha

    TPZMaterial * BCond1 = material->CreateBC(material, fmatBCright, fneumann, val1, val2); //Cria material que implementa a condicao de contorno superior
    BCond1->SetForcingFunction(Loads::Sol_exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond1); //Insere material na malha

    TPZMaterial * BCond2 = material->CreateBC(material, fmatBCtop, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno esquerda
    BCond2->SetForcingFunction(Loads::Sol_exact,fIntOrder);
    cmesh->InsertMaterialObject(BCond2); //Insere material na malha

    TPZMaterial * BCond3 = material->CreateBC(material, fmatBCleft, fdirichlet, val1, val2); //Cria material que implementa a condicao de contorno direita
    BCond3->SetForcingFunction(Loads::Sol_exact,fIntOrder);
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