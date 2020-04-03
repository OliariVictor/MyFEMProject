//
// Created by victor on 19/02/2020.
//

#include "MixedHdiv.h"

DarcyMixedHdiv::DarcyMixedHdiv(int pOrder): fnx(2,1), fhx(2,1){
    fdim = 2;
    fPOrder = pOrder;
    fdebug = false;
}

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

void DarcyMixedHdiv::ErrorRate(int numRefinement,TLaplaceExample1 &Laplace,TPZFMatrix<REAL>  K, TPZFMatrix<REAL>  invK, bool postproc){
    if (numRefinement < 0 || numRefinement > 5) {std::cout << "Invalid refinement number"; DebugStop();}

    TPZVec<REAL> x0(3,0),x1(3,0);
    for(int i =0; i <fdim;i++)
        x1[i] = fhx[i];

    Geometry geom(fdim);

    int numErrors = 5;
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)

    //Clearing file
    if( remove( "Erro.txt" ) != 0 || remove( "ErroHybrid.txt" ) != 0 )
        perror( "Error deleting file" );
    else
        puts( "Error log successfully deleted" );

    //Defining outstream file
    ofstream ErroOut("Erro.txt",std::ofstream::app);
    ofstream ErroOutHybrid("ErroHybrid.txt",std::ofstream::app);
    ErroOut << "----------COMPUTED ERRORS----------\n"; ErroOut.flush();
    ErroOutHybrid << "----------COMPUTED ERRORS (HYBRID)----------\n"; ErroOutHybrid.flush();

    TPZVec<REAL> ErrorsLog(numErrors,-1), ErrorsLogHybrid(numErrors,-1),rate(3,-1);
    REAL hLog = -1, h;

    bool last = false;
    int exp = 1; //Initial exponent of mesh refinement (numElem = 2*2^exp)
    for(int i =0 ; i <numRefinement+1 ; i++){
        fnx[0] = fnx[1] = exp;
        h = 1./exp;

        TPZGeoMesh* gmesh = geom.gmeshGenGrid(fnx,x0,x1,1);

        Computational comp;
        TPZCompMesh *cmesh_f = comp.CMesh_flux(gmesh,fPOrder);
        TPZCompMesh *cmesh_p = comp.CMesh_p(gmesh,fPOrder);
        TPZMultiphysicsCompMesh *cmesh_m = comp.CMesh_m(gmesh,fPOrder,Laplace,K,invK);

        TPZManVector<int> active(2, 1);
        TPZManVector<TPZCompMesh *, 2> meshvector(2);
        meshvector[0] = cmesh_f;
        meshvector[1] = cmesh_p;
        cmesh_m->BuildMultiphysicsSpace(active,meshvector);
        cmesh_m->LoadReferences();
        cmesh_m->InitializeBlock();

        std::ofstream oGmesh("gmeshMixed.txt");
        gmesh->Print(oGmesh);

        std::ofstream pressure("cmesh_p.txt");
        cmesh_p->Print(pressure);

        std::ofstream flux("cmesh_f.txt");
        cmesh_f->Print(flux);

        if(i == numRefinement && postproc) last = true;

        if(1) {
            SolveMixedProblem(cmesh_m, optimizeBandwidth, Laplace, ErrorsLog, rate, gmesh, hLog, h,last);
            std::ofstream mixed("MultiphysicsCompMesh.txt");
            cmesh_m->ShortPrint(mixed);
            std::ofstream mixedP("MultiphysicsCompMeshPrint.txt");
            cmesh_m->Print(mixedP);
        }

        if(1){
            TPZHybridizeHDiv hybrid;
            TPZMultiphysicsCompMesh *HybridMesh = hybrid.Hybridize(cmesh_m);
            HybridMesh->CleanUpUnconnectedNodes();
            HybridMesh->AdjustBoundaryElements();

            std::ofstream oGmeshH("gmeshHybrid.txt");
            gmesh->Print(oGmeshH);
            ofstream gvtk("gm.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, gvtk);

            SolveHybridProblem(HybridMesh, optimizeBandwidth, Laplace, ErrorsLogHybrid, rate, gmesh, hLog, h,last);

            std::ofstream hybred("HybridMultiphysicsCompMesh.txt");
            HybridMesh->ShortPrint(hybred);
            std::ofstream hybredP("HybridMultiphysicsCompMeshPrint.txt");
            HybridMesh->Print(hybredP);
        }
        hLog = h;
        exp *=2;
    }
}

void DarcyMixedHdiv::SolveHybridProblem(TPZCompMesh *mesh, bool optBW, TLaplaceExample1 &Laplace, TPZVec<REAL> &errorsLog, TPZVec<REAL> &rt, TPZGeoMesh* gmesh, REAL &hLog, REAL &h,bool &last){
    ofstream erroOutput("ErroHybrid.txt",std::ofstream::app);

    TPZAnalysis an(mesh, optBW); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineStructMatrix matskl(mesh); //caso simetrico ***

    int numthreads = 0;
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    an.Solve();

    //Calculo do erro
    std::cout << "Computing Error " << std::endl;
    TPZManVector<REAL,6> Errors;

    an.SetExact(Laplace.ExactSolution());
    bool store_errors = false;
    an.PostProcessError(Errors, store_errors, erroOutput);

    if(errorsLog[0] !=-1){
        for(int j =0; j < 3 ; j++) {
            rt[j] = (log10(Errors[j]) - log10(errorsLog[j]))/(log10(h)-log10(hLog));
            erroOutput << "rate " << j << ": " << rt[j] << std::endl;
        }
    }

    erroOutput << "h = " << 1./(fnx[0]) << std::endl;

    if(last) {
        int dim = gmesh->Dimension();
        string Name = "Mixed.vtk";
        PostProcess(Name,dim,an);
    }

    std::cout << "FINISHED!" << std::endl;

    for(int i =0; i < Errors.size();i++) errorsLog[i] = Errors[i];
    Errors.clear();
}

void DarcyMixedHdiv::SolveMixedProblem(TPZCompMesh *mesh, bool optBW, TLaplaceExample1 &Laplace, TPZVec<REAL> &errorsLog, TPZVec<REAL> &rt, TPZGeoMesh* gmesh, REAL &hLog, REAL &h,bool &last){
    ofstream erroOutput("Erro.txt",std::ofstream::app);

    TPZAnalysis an(mesh, optBW); //Cria objeto de análise que gerenciará a analise do problema
    TPZSkylineStructMatrix matskl(mesh); //caso simetrico ***

    int numthreads = 0;
    matskl.SetNumThreads(numthreads);
    an.SetStructuralMatrix(matskl);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
    an.Solve();

    //Calculo do erro
    std::cout << "Computing Error " << std::endl;
    TPZManVector<REAL,6> Errors;

    an.SetExact(Laplace.ExactSolution());
    bool store_errors = false;
    an.PostProcessError(Errors, store_errors, erroOutput);

    if(errorsLog[0] !=-1){
        for(int j =0; j < 3 ; j++) {
            rt[j] = (log10(Errors[j]) - log10(errorsLog[j]))/(log10(h)-log10(hLog));
            erroOutput << "rate " << j << ": " << rt[j] << std::endl;
        }
    }

    erroOutput << "h = " << 1./(fnx[0]) << std::endl;

    if(last) {
        int dim = gmesh->Dimension();
        string Name = "Hybrid.vtk";
        PostProcess(Name,dim,an);
    }

    std::cout << "FINISHED!" << std::endl;

    for(int i =0; i < Errors.size();i++) errorsLog[i] = Errors[i];
    Errors.clear();
}

void DarcyMixedHdiv::PostProcess(string fileName, int dim, TPZAnalysis &analysis){

    std::cout << "Post Processing " << std::endl;
    std::string plotfile(fileName);
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("ExactPressure");
    vecnames.Push("Flux");
    vecnames.Push("ExactFlux");

    int postProcessResolution = 4; //  keep low as possible

    analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    analysis.PostProcess(postProcessResolution,dim);
}

void DarcyMixedHdiv::ErrorRateProperFunc(int numRefinement,void (*f_source)(const TPZVec<REAL> &x, TPZVec<STATE> &val),void (*Sol_Exact)(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf)){
    if (numRefinement < 0 || numRefinement > 5) {std::cout << "Invalid refinement number"; DebugStop();}

    TPZVec<REAL> x0(3,0),x1(3,0);
    for(int i =0; i <fdim;i++)
        x1[i] = fhx[i];

    Geometry geom(fdim);

    int numthreads = 0;
    bool optimizeBandwidth = true; //Impede a renumeração das equacoes do problema (para obter o mesmo resultado do Oden)

    //Clearing file
    if( remove( "Erro.txt" ) != 0 )
        perror( "Error deleting file" );
    else
        puts( "Error log successfully deleted" );

    //Defining outstream file
    ofstream ErroOut("Erro.txt",std::ofstream::app);
    ErroOut << "----------COMPUTED ERRORS----------\n";

    int exp = 1; //Initial exponent of mesh refinement (numElem = 2*2^exp)
    for(int i =0 ; i <numRefinement+1 ; i++){
        fnx[0] = fnx[1] = exp;

        TPZGeoMesh* gmesh = geom.gmeshGenGrid(fnx,x0,x1,1);

        Computational comp;
        TPZCompMesh *cmesh_f = comp.CMesh_flux(gmesh,fPOrder);
        TPZCompMesh *cmesh_p = comp.CMesh_p(gmesh,fPOrder);
        TPZCompMesh *cmesh_m = comp.CMeshProperFunc_m(gmesh,fPOrder,f_source,Sol_Exact);

        TPZManVector<TPZCompMesh *, 2> meshvector(2);
        meshvector[0] = cmesh_f;
        meshvector[1] = cmesh_p;
        TPZBuildMultiphysicsMesh::AddElements(meshvector, cmesh_m);
        TPZBuildMultiphysicsMesh::AddConnects(meshvector, cmesh_m);

        //Resolvendo o Sistema:
        TPZAnalysis an(cmesh_m, optimizeBandwidth); //Cria objeto de análise que gerenciará a analise do problema
        TPZSkylineStructMatrix matskl(cmesh_m); //caso simetrico ***

        matskl.SetNumThreads(numthreads);
        an.SetStructuralMatrix(matskl);
        TPZStepSolver<STATE> step;
        step.SetDirect(ELDLt);
        an.SetSolver(step);
        an.Assemble();//Assembla a matriz de rigidez (e o vetor de carga) global
        an.Solve();

        //Calculo do erro
        std::cout << "Comuting Error " << std::endl;
        TPZManVector<REAL,3> Errors;

        an.SetExact(Sol_Exact);
        bool store_errors = false;
        an.PostProcessError(Errors, store_errors, ErroOut);
        ErroOut << "h = " << 1./(fnx[0]) << std::endl;
        std::cout << Errors[0] << "; " << Errors[1] << "; " << Errors[2] << std::endl;
        Errors.clear();

        exp *=2;
    }
}

void DarcyMixedHdiv::SolveMixedHdiv(TLaplaceExample1 &Laplace){
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
    TPZCompMesh *cmesh_m = comp.CMesh_m(gmesh,fPOrder,Laplace);

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
    cmesh_m->ShortPrint(filecmfp);
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
    an.SetExact(Laplace.ExactSolution());
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

void DarcyMixedHdiv::SolveMixedHdivProperFunc(void (*f_source)(const TPZVec<REAL> &x, TPZVec<STATE> &val),void (*Sol_Exact)(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf)){
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
    TPZCompMesh *cmesh_m = comp.CMeshProperFunc_m(gmesh,fPOrder,f_source,Sol_Exact);

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
    an.SetExact(Sol_Exact);
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