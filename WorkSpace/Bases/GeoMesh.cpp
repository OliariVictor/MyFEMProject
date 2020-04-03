//
// Created by victor on 18/02/2020.
//

#include "GeoMesh.h"

Geometry::Geometry(){
    fdim = 2;
    fmatId = 1;

    fmatBCbottom = -1;
    fmatBCright = -2;
    fmatBCtop = -3;
    fmatBCleft = -4;
}

Geometry::Geometry(int dim){
    fdim = dim;
    fmatId = 1;

    fmatBCbottom = -1;
    fmatBCright = -2;
    fmatBCtop = -3;
    fmatBCleft = -4;
}

TPZGeoMesh* Geometry::gmeshHardWay(int nx, int ny, double hx, double hy){
    int i,j;
    int fmatID = 1;
    int64_t id, index;


    //Criando malha geométrica, nós e elementos.
    //Inserindo nós e elementos no objeto malha:

    TPZGeoMesh *gmesh = new TPZGeoMesh();
    gmesh->SetDimension(fdim);

    //Vetor auxiliar para armazenar coordenadas:

    TPZVec <REAL> coord (3,0.);


    //Inicialização dos nós:

    for(i = 0; i < ny; i++){
        for(j = 0; j < nx; j++){
            id = i*nx + j;
            coord[0] = (j)*hx/(nx - 1);
            coord[1] = (i)*hy/(ny - 1);
            //using the same coordinate x for z
            coord[2] = 0.;
            //cout << coord << endl;
            //Get the index in the mesh nodes vector for the new node
            index = gmesh->NodeVec().AllocateNewElement();
            //Set the value of the node in the mesh nodes vector
            gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
        }
    }

    TPZVec <int64_t> connectD(3,0);
    TPZVec <int64_t> connectU(3,0);


    //Conectividade dos elementos:

    for(i = 0; i < (ny - 1); i++){
        for(j = 0; j < (nx - 1); j++){
            index = (i)*(nx - 1)+ (j);
            connectD[0] = (i)*ny + (j);
            connectD[1] = connectD[0]+1;
            connectD[2] = connectD[1]+nx;
            gmesh->CreateGeoElement(ETriangle,connectD,fmatId,id);

            connectU[0] = connectD[0];
            connectU[1] = connectD[1]+nx;
            connectU[2] = connectD[0]+nx;
            gmesh->CreateGeoElement(ETriangle,connectU,fmatId,id);

            id++;
        }
    }
    gmesh->BuildConnectivity();

    //Finding each element boundary
    //1. Iterate through each element;
    //2. On each element determine weather a node belongs to an internal boundary or to the skeleton;
    //3. Create a vector with the bottom, right, top and left nodes;
    //4. Create a Boundary Element with the corresponding nodes

    int64_t numEl = gmesh->NElements();
    for(int64_t iel = 0; iel < numEl ;iel++){
        TPZGeoEl *geoel = gmesh->Element(iel);
        int64_t nnodes = geoel->NNodes();
        TPZVec<int> nodeInd(nnodes,-1);

        for(int ni = 0; ni < nnodes; ni++) nodeInd[ni] = geoel->NodeIndex(ni);

        TPZManVector<REAL,2> bottom(0); bottom.resize(0);
        TPZManVector<REAL,2> right(0); right.resize(0);
        TPZManVector<REAL,2> top(0); top.resize(0);
        TPZManVector<REAL,2> left(0); left.resize(0);

        for(int64_t in = 0; in < nnodes ; in++){
            TPZManVector<REAL,3> coord(3,0);
            gmesh->NodeVec()[nodeInd[in]].GetCoordinates(coord);

            if(coord[1] == 0){ //bottomside
                bottom.resize(bottom.size()+1);
                bottom[bottom.size()-1] = in;
            }if(coord[1] == (ny-1)*hy){//topside
                top.resize(top.size()+1);
                top[top.size()-1] = in;
            }
            if(coord[0] == 0){//leftside
                left.resize(left.size()+1);
                left[left.size()-1] = in;
            }
            if(coord[0] == hx*(nx-1)){//rightside
                right.resize(right.size()+1);
                right[right.size()-1] = in;
            }
        }

        TPZVec<int64_t> nodeIds(2,-1);
        if(bottom.size() == 2){
            for(int i =0; i <2;i++) nodeIds[i] = nodeInd[bottom[i]];
            int sidenum = geoel->WhichSide(nodeIds);
            TPZGeoElSide geoSide(geoel,sidenum);
            TPZGeoElBC geoBC(geoSide,-1);
            //gmesh->CreateGeoElement(EOned, nodeIds,-1,id);
        }
        if(top.size() == 2){
            for(int i =0; i <2;i++) nodeIds[i] = nodeInd[top[i]];
            int sidenum = geoel->WhichSide(nodeIds);
            TPZGeoElSide geoSide(geoel,sidenum);
            TPZGeoElBC geoBC(geoSide,-2);
            //gmesh->CreateGeoElement(EOned, nodeIds,-1,id);
        }
        if(left.size() == 2){
            for(int i =0; i <2;i++) nodeIds[i] = nodeInd[left[i]];
            int sidenum = geoel->WhichSide(nodeIds);
            TPZGeoElSide geoSide(geoel,sidenum);
            TPZGeoElBC geoBC(geoSide,-3);
            //gmesh->CreateGeoElement(EOned, nodeIds,-1,id);
        }
        if(right.size() == 2){
            for(int i =0; i <2;i++) nodeIds[i] = nodeInd[right[i]];
            int sidenum = geoel->WhichSide(nodeIds);
            TPZGeoElSide geoSide(geoel,sidenum);
            TPZGeoElBC geoBC(geoSide,-4);
            //gmesh->CreateGeoElement(EOned, nodeIds,-1,id);
        }

        bottom.resize(0);top.resize(0);right.resize(0);left.resize(0);
    }
    TPZCheckGeom check(gmesh);
    check.CheckUniqueId();

    gmesh->BuildConnectivity();

    //Impressão da malha geométrica:
    return gmesh;
}

TPZGeoMesh* Geometry::gmeshGenGrid(TPZVec<int> &nx, TPZVec<REAL> &x0, TPZVec<REAL> &x1,int matId){
    TPZGenGrid2D genGrid(nx,x0,x1);
    MMeshType shape = MMeshType::EQuadrilateral;
    genGrid.SetElementType(shape);

    //TPZGenGrid genGrid(nx,x0,x1);
    //genGrid.SetElementType(ETriangle);

    TPZGeoMesh *gmesh = new TPZGeoMesh(); gmesh->SetDimension(2);
    genGrid.Read(gmesh,matId);

    genGrid.SetBC(gmesh,4,fmatBCbottom);
    genGrid.SetBC(gmesh,5,fmatBCright);
    genGrid.SetBC(gmesh,6,fmatBCtop);
    genGrid.SetBC(gmesh,7,fmatBCleft);

    gmesh->BuildConnectivity();
    return gmesh;
}

void Geometry::TestHWvsGG(){
    TPZVec<int> nx(2,0), hx(2,0);
    nx[0] = 2; nx[1] = 2;
    hx[0] = 1; hx[1] = 1;

    TPZGeoMesh* gelHWzero = gmeshHardWay(nx[0],nx[1], hx[0], hx[1]);
    TPZGeoMesh* gelHWone = gmeshHardWay(nx[0]+1,nx[1]+1, hx[0], hx[1]);
    TPZGeoMesh* gelHWtwo = gmeshHardWay(nx[0]+2,nx[1]+2, hx[0], hx[1]);

    TPZVec<REAL> x0 ={0,0,0},x1 = {hx[0],hx[1],0};
    TPZGeoMesh* gelGGzero = gmeshGenGrid(nx,x0,x1,1); nx[0]++;nx[1]++; std::cout << nx[0];
    TPZGeoMesh* gelGGone = gmeshGenGrid(nx,x0,x1,1);nx[0]++;nx[1]++;
    TPZGeoMesh* gelGGtwo = gmeshGenGrid(nx,x0,x1,1);

    ofstream HW0("gelHWzero.vtk"),HW1("gelHWone.vtk"),HW2("gelHWtwo.vtk");
    ofstream GG0("gelGGzero.vtk"),GG1("gelGGone.vtk"),GG2("gelGGtwo.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gelHWzero, HW0);TPZVTKGeoMesh::PrintGMeshVTK(gelHWone, HW1);TPZVTKGeoMesh::PrintGMeshVTK(gelHWtwo, HW2);
    TPZVTKGeoMesh::PrintGMeshVTK(gelGGzero, GG0);TPZVTKGeoMesh::PrintGMeshVTK(gelGGone, GG1);TPZVTKGeoMesh::PrintGMeshVTK(gelGGtwo, GG2);
}