#include "simplicialcomplex.h"

Vertex SimplicialComplex::barycenter(const explicitS &s){

    vector<int> vs = getTopSimplex(s).getVertices();

    vector<float> f(3,0);

    for(auto v : vs){
        for(int i=0; i<3; i++){
            f[i]+=getVertex(v).getCoordinate(i);
        }
    }

    for(int i=0; i<3; i++)
        f[i] = f[i]/vs.size();

    return Vertex(f);
}

Vertex SimplicialComplex::barycenter(const implicitS &s){

    vector<int> vs = s.getConstVertices();

    vector<float> f(3,0);

    for(auto v : vs){
        for(int i=0; i<3; i++){
            f[i]+=getVertex(v).getCoordinate(i);
        }
    }

    for(int i=0; i<3; i++)
        f[i] = f[i]/vs.size();

    return Vertex(f);
}


vector<float> SimplicialComplex::vectorSum(const vector<float>& v1, const vector<float>& v2){

    assert(v1.size()==3);
    assert(v2.size()==3);

    vector<float> retV={0,0,0};
    for(unsigned int i=0;i<v1.size();i++)
        retV[i]=(v1[i]+v2[i]);

    return retV;
}

vector<float> SimplicialComplex::vectorSubtr(const vector<float>& v1, const vector<float>& v2){

    assert(v1.size()==3);
    assert(v2.size()==3);

    vector<float> retV={0,0,0};
    for(unsigned int i=0;i<v1.size();i++)
        retV[i]=(v1[i]-v2[i]);

    return retV;
}

vector<float> SimplicialComplex::scalarProd(const vector<float> &v1, const vector<float> &v2){

    assert(v1.size()==3);
    assert(v2.size()==3);

    vector<float> retV={0,0,0};

    retV[0] = (v1[1]*v2[2]-v1[2]*v2[1]);
    retV[1] = (v1[2]*v2[0]-v1[0]*v2[2]);
    retV[2] = (v1[0]*v2[1]-v1[1]*v2[0]);

    return retV;
}

vector<float> SimplicialComplex::prod(const vector<float> &v1, float prod){
    assert(v1.size()==3);
    vector<float> retV={0,0,0};
    for(int i=0; i<3; i++)
        retV[i]=v1[i]*prod;
    return retV;
}

float SimplicialComplex::dotProd(const vector<float> &v1, const vector<float> &v2){
     assert(v1.size()==3);
     assert(v2.size()==3);
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

float SimplicialComplex::dist(const vector<float> &v1, const vector<float> &v2){
    assert(v1.size()==3);
    assert(v2.size()==3);
    return sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2)+pow(v1[2]-v2[2],2));
}

float SimplicialComplex::norm(const vector<float> &v1){
    assert(v1.size()==3);
    return sqrt(pow(v1[0],2)+pow(v1[1],2)+pow(v1[2],2));
}

