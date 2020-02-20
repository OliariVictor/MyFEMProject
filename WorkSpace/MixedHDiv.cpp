//
// Created by victor on 19/02/2020.
//

#include "MixedHdiv.h"

DarcyMixedHdiv::DarcyMixedHdiv(TPZVec<int> nx,TPZVec<int> hx,int dim,int pOrder, bool debug){
    fdim = dim;
    if(nx.size() != dim || hx.size() != dim){
        std::cout<< "MixedHDiv::DarcyMixedHdiv(TPZVec<int> nx,TPZVec<int> hx,int dim,int pOrder) \n\tsize of nx/hx not compatible with dimension";
        DebugStop();
    }
    fnx = nx;
    fhx = hx;
    fPOrder = pOrder;
    fdebug = debug;
}

DarcyMixedHdiv::~DarcyMixedHdiv(){

}

void DarcyMixedHdiv::SolveMixedHdiv(){
    //Definition of x
    TPZVec<REAL> x0(3,0),x1(3,0);
    for(int i =0; i <fdim;i++)
        x1[i] = fhx[i];

    if(fdebug){
        std::cout << "fdebug = true\n";
        #ifndef DEBUG
        #define DEBUG
        #endif
    }
    //Creation of geometric mesh
    Geometry geom(fdim);
    TPZGeoMesh* gmesh = geom.gmeshGenGrid(fnx,x0,x1,1);

#ifdef DEBUG
    //Printing Geometrical mesh
    ofstream gmName("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, gmName);
#endif

    Computational comp;
    TPZCompMesh *cmesh_f = comp.CMesh_flux(gmesh,fPOrder);
    TPZCompMesh *cmesh_p = comp.CMesh_p(gmesh,fPOrder);
    TPZCompMesh *cmesh_m = comp.CMesh_m(gmesh,fPOrder);

#ifdef DEBUG
    //Printing Computational mesh
    std::ofstream filecf("MalhaC_v.txt"); //Impressão da malha computacional da velocidade (formato txt)
    std::ofstream filecp("MalhaC_p.txt"); //Impressão da malha computacional da pressão (formato txt)
    std::ofstream filecm("MalhaC_m.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_f->Print(filecf);
    cmesh_p->Print(filecp);
    cmesh_m->Print(filecm);
#endif

    TPZManVector<TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmesh_f;
    meshvector[1] = cmesh_p;
    TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);

#ifdef DEBUG
    //Printing Multiphisics mesh after the insertion of the atomic meshes
    std::ofstream filecmfp("MalhaC_m_add_fp.txt"); //Impressão da malha computacional multifísica (formato txt)
    cmesh_m->Print(filecmfp);
#endif

    //Resolvendo o Sistema:
    int numthreads = 0;
    bool optimizeBandwidth = false; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)

    TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineStructMatrix matskl(cmesh_m); //caso simetrico ***

    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);


    std::cout << "Assemble matrix with NDoF = " << cmesh_m->NEquations() << std::endl;

    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global

#ifdef DEBUG
    //Banded Global matrix
    TPZFMatrix<REAL> fillin;
    int resolution = 100;
    cmesh_m->ComputeFillIn(resolution, fillin);
    std::string out("matrix_native.vtk");
    VisualMatrix(fillin, out);
#endif
    std::cout << "Solving Matrix " << std::endl;

    an.Solve();

    //Calculo do erro
    std::cout << "Comuting Error " << std::endl;
    TPZManVector<REAL,3> Errors;
    ofstream ErroOut("Erro.txt");
    an.SetExact(Loads::Sol_exact);
    bool store_errors = false;
    an.PostProcessError(Errors, store_errors, ErroOut);



    //Pós-processamento (paraview):
    std::cout << "Post Processing " << std::endl;
    std::string plotfile("DarcyP.vtk");
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");
    //        vecnames.Push("V_exactBC");


    int postProcessResolution = 4; //  keep low as possible

    int dim = gmesh->Dimension();
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(postProcessResolution,dim);

    std::cout << "FINISHED!" << std::endl;
}

