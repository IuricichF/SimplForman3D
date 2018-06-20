#include "gradientencoding.h"

GradientEncoding::GradientEncoding(SimplicialComplex& sc){

    edgeLocalFrames = vector<bool>(2*sc.getTopSimplexesNum(1),false);
    triangleLocalFrames = vector<unsigned short>(sc.getTopSimplexesNum(2),0);
    tetraLocalFrames = vector<unsigned short>(sc.getTopSimplexesNum(3),0);

    FILE* lookup;

    //read the valid configurations from file (2D case)
    lookup = fopen("table2D.txt","r");

    if(lookup != NULL){
        unsigned int arrow,triCase;
        triVect = vector<TriangleGradient>(42);
        for(int i=0; i<42; i++){
            fscanf(lookup,"%u %u\n", &arrow, &triCase);

            TriangleGradient tg = TriangleGradient(arrow);
            triLookup[tg]=i;
            triVect[i]=tg;
        }
        fclose(lookup);
    }
    else
        cout << "Error - file ../data/table2D.txt cannot be found" << endl;

    //read the valid configurations from file (3D case)
    lookup = fopen("table3D.txt","r");

    if(lookup != NULL){
        unsigned int arrow,tetCase;
        tetraVect = vector<TetraGradient>(51030);
        for(int i=0; i<51030; i++){
            fscanf(lookup,"%u %u\n", &arrow, &tetCase);

            TetraGradient tg = TetraGradient(arrow);
            tetraLookup[tg]=i;
            tetraVect[i]=tg;
        }
        fclose(lookup);
    }
    else
        cout << "Error - file ../data/table3D.txt cannot be found" << endl;
}

TriangleGradient& GradientEncoding::expandTriangleGrad(unsigned short grad){
    return triVect[grad];
}

TetraGradient& GradientEncoding::expandTetraGrad(unsigned short grad){
    return tetraVect[grad];
}

void GradientEncoding::compressGradTri(TriangleGradient const& tri, explicitS const& simpl){
    triangleLocalFrames[simpl.getIndex()]=triLookup[tri];
}

void GradientEncoding::compressGradTet(TetraGradient const& tet, explicitS const& simpl){
    tetraLocalFrames[simpl.getIndex()]=tetraLookup[tet];
}

SBitSet GradientEncoding::getBitSet(const implicitS& s1,const implicitS& s2){

    bitset<4> s;

    for(int i=0; i<s2.getConstVertices().size(); i++){
        for(auto v : s1.getConstVertices()){
            if(s2.getConstVertices()[i] == v){
                s[i]=1;
            }
        }
    }

    return SBitSet(s.to_ulong());
}

implicitS GradientEncoding::bitSetToImplicitS(const SBitSet& s1,const implicitS& s2){

    bitset<4> sbit(s1);
    vector<int> vertices;
    for(int i=0; i<4; i++){
        if(sbit[i]){
            vertices.push_back(s2.getConstVertices()[i]);
        }
    }

    return implicitS(vertices);
}

bool GradientEncoding::isPaired(const implicitS& simpl, explicitS& top, SimplicialComplex& sc){


    uint topIndex = top.getIndex();
    uint topDim = top.getDim();

    vector<int> vsimpl = simpl.getConstVertices();
    vector<int> vtop = sc.getTopSimplex(top).getVertices();

    SBitSet vbit = getBitSet(simpl,sc.toImplicit(top));

    switch(topDim){

        case 1:
            if(topDim == 1){
                return edgeLocalFrames[topIndex] || edgeLocalFrames[topIndex+1];
            }

            if(vsimpl[0] == vtop[0])
                return edgeLocalFrames[topIndex];

            if(vsimpl[0] == vtop[1])
                return edgeLocalFrames[topIndex+1];
            else{
                cout << "ERROR IN ISPAIRED1" << endl;
                return false;
            }
            break;

        case 2:

            switch(simpl.getDim()){
                case 0:
                    return !expandTriangleGrad(triangleLocalFrames[topIndex]).is_vertex_unpaired(vbit);
                case 1:
                    return !expandTriangleGrad(triangleLocalFrames[topIndex]).is_edge_unpaired(vbit);
                case 2:
                    return !expandTriangleGrad(triangleLocalFrames[topIndex]).is_triangle_unpaired();
                default:
                    cout << "ERROR IN ISPAIRED2" << endl;
                    break;
            }
            break;
        case 3:
            switch(simpl.getDim()){
                case 0:
                    return !expandTetraGrad(tetraLocalFrames[topIndex]).is_vertex_unpaired(vbit);
                case 1:
                    return !expandTetraGrad(tetraLocalFrames[topIndex]).is_edge_unpaired(vbit);
                case 2:
                    return !expandTetraGrad(tetraLocalFrames[topIndex]).is_face_unpaired(vbit);
                case 3:
                    return !expandTetraGrad(tetraLocalFrames[topIndex]).is_tetra_unpaired();
                default:
                    cout << "ERROR IN ISPAIRED3" << endl;
                    break;
        }
        break;
        default:
            cout << " ERROR IN ISPAIRED dimension" << endl;
    }

    return false;
}


void GradientEncoding::pair(const implicitS& big,const implicitS& small, explicitS& top, SimplicialComplex& sc){

    assert(small.getDim() < big.getDim());

    uint smallDim = small.getDim();
    SBitSet smallBit = getBitSet(small,sc.toImplicit(top));

    if(top.getDim() == 1){

        switch (smallBit){
            case 1:
                edgeLocalFrames[top.getIndex()]=1;
                break;
            case 2:
                edgeLocalFrames[top.getIndex()+1]=1;
                break;
            default:
                cout << " ERROR IN PAIR 1" << endl;
                break;
        }
    }

    else if(top.getDim() == 2){

        TriangleGradient tg = expandTriangleGrad(triangleLocalFrames[top.getIndex()]);
        switch(smallDim){
            case 0:
                tg.setVE(getBitSet(small,sc.toImplicit(top)),getBitSet(big,sc.toImplicit(top)));
                break;
            case 1:
                tg.setEF(getBitSet(small,sc.toImplicit(top)));
                break;
            default:
            cout << "ERROR IN PAIR 2" << endl;
            break;
        }
        compressGradTri(tg,top);
    }

    else{

        TetraGradient tg = expandTetraGrad(tetraLocalFrames[top.getIndex()]);
        switch(smallDim){
            case 0:
                tg.setVE(getBitSet(small,sc.toImplicit(top)),getBitSet(big,sc.toImplicit(top)));
                break;
            case 1:
                tg.setEF(getBitSet(small,sc.toImplicit(top)),getBitSet(big,sc.toImplicit(top)));
                break;
            case 2:
                tg.setFT(getBitSet(small,sc.toImplicit(top)));
                break;
            default:
                cout << "ERROR IN PAIR 3" << endl;
                break;
        }
        compressGradTet(tg,top);
    }

}


void GradientEncoding::free(const implicitS& big,const implicitS& small, explicitS& top, SimplicialComplex& sc){

    assert(small.getDim() < big.getDim());
    uint smallDim = small.getDim();

    if(top.getDim() == 1){

        switch (getBitSet(small,big)){
            case 1:
                edgeLocalFrames[top.getIndex()]=0;
                break;
            case 2:
                edgeLocalFrames[top.getIndex()+1]=0;
                break;
            default:
                cout << " ERROR IN PAIR 1" << endl;
                break;
        }
    }

    else if(top.getDim() == 2){

        TriangleGradient tg = expandTriangleGrad(triangleLocalFrames[top.getIndex()]);
        switch(smallDim){
            case 0:
                tg.freeVE(getBitSet(small,sc.toImplicit(top)),getBitSet(big,sc.toImplicit(top)));
                break;
            case 1:
                tg.freeEF(getBitSet(small,sc.toImplicit(top)));
                break;
            default:
            cout << "ERROR IN PAIR 2" << endl;
            break;
        }
        compressGradTri(tg,top);
    }

    else{

        TetraGradient tg = expandTetraGrad(tetraLocalFrames[top.getIndex()]);
        switch(smallDim){
            case 0:
                tg.freeVE(getBitSet(small,sc.toImplicit(top)),getBitSet(big,sc.toImplicit(top)));
                break;
            case 1:
                tg.freeEF(getBitSet(small,sc.toImplicit(top)),getBitSet(big,sc.toImplicit(top)));
                break;
            case 2:
                tg.freeFT(getBitSet(small,sc.toImplicit(top)));
                break;
            default:
                cout << "ERROR IN PAIR 3" << endl;
                break;
        }
        compressGradTet(tg,top);
    }

}


bool GradientEncoding::getPair(const implicitS& simpl, implicitS& next, explicitS& top, SimplicialComplex& sc){

    uint topIndex = top.getIndex();

    SBitSet pair=EMPTY;
    switch(top.getDim()){

        case 1:
            switch(getBitSet(simpl,sc.toImplicit(top))){
                case 1:
                    if(edgeLocalFrames[topIndex]){
                        next = sc.toImplicit(top); // return edge
                    }
                    return edgeLocalFrames[topIndex];
                case 2:
                    if(edgeLocalFrames[topIndex+1]){
                        next = sc.toImplicit(top); // return edge
                    }
                    return edgeLocalFrames[topIndex+1];
                case 3:
                    if(edgeLocalFrames[topIndex]){
                        vector<int> s;
                        s.push_back(simpl.getConstVertices()[0]);
                        next = implicitS(s);
                    }
                    else if(edgeLocalFrames[topIndex+1]){
                        vector<int> s;
                        s.push_back(simpl.getConstVertices()[1]);
                        next = implicitS(s);
                    }

                    return edgeLocalFrames[topIndex] || edgeLocalFrames[topIndex+1];

            }
            break;
        case 2:
            switch(simpl.getDim()){
                case 0:
                    pair = expandTriangleGrad(triangleLocalFrames[topIndex]).get_vertex_pair(getBitSet(simpl,sc.toImplicit(top)));
                    break;
                case 1:
                    pair = expandTriangleGrad(triangleLocalFrames[topIndex]).get_edge_pair(getBitSet(simpl,sc.toImplicit(top)));
                    break;
                case 2:
                    pair = expandTriangleGrad(triangleLocalFrames[topIndex]).get_triangle_pair();
                    break;
                default:
                    cout << "ERROR IN ISPAIRED2" << endl;
                    break;
            }
            next = bitSetToImplicitS(pair,sc.toImplicit(top));
            return pair != EMPTY;
            break;
        case 3:

            switch(simpl.getDim()){
                case 0:
                    pair = expandTetraGrad(tetraLocalFrames[topIndex]).get_vertex_pair(getBitSet(simpl,sc.toImplicit(top)));
                    break;
                case 1:
                    pair = expandTetraGrad(tetraLocalFrames[topIndex]).get_edge_pair(getBitSet(simpl,sc.toImplicit(top)));
                    break;
                case 2:
                    pair = expandTetraGrad(tetraLocalFrames[topIndex]).get_face_pair(getBitSet(simpl,sc.toImplicit(top)));
                    break;
                case 3:
                    pair = expandTetraGrad(tetraLocalFrames[topIndex]).get_tetra_pair();
                    break;
                default:
                    cout << "ERROR IN ISPAIRED3" << endl;
                    break;
            }
            next = bitSetToImplicitS(pair,sc.toImplicit(top));
            return pair != EMPTY;
            break;
            default:
                cout << " ERROR IN ISPAIRED dimension" << endl;
                break;
    }

    return false;
}
