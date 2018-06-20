#include "topsimplex.h"

//--------------------------------------
//implicitS implementation

implicitS::implicitS(){vertices=vector<int>();}

implicitS::implicitS(int v) {
    vertices = vector<int>(1);
    vertices[0]=v;
}

implicitS::implicitS(const vector<int>& simplexV) : vertices(simplexV) {

}

uint implicitS::getDim() const {assert(vertices.size()>0); return vertices.size()-1;}

vector<int>* implicitS::getVertices(){ return &vertices;}

vector<int> const& implicitS::getConstVertices() const{ return vertices;}

bool implicitS::operator!=(const implicitS& simpl) const{
    return !(*this==simpl);
}

bool implicitS::contains(const implicitS &simpl){

    vector<int> sVerts = simpl.getConstVertices();
    if(vertices.size() <= sVerts.size())
        return false;

    for(auto v : sVerts){
        if(find(vertices.begin(), vertices.end(), v) == vertices.end()){
            return false;
        }
    }

    return true;
}

bool implicitS::operator==(const implicitS& simpl) const{
    if(this->getDim() == simpl.getDim())
        return equal(this->getConstVertices().begin(), this->getConstVertices().end(),simpl.getConstVertices().begin());
    return false;
}

bool implicitS::operator <(const implicitS& simpl) const{
    if(getDim() == simpl.getDim()){

        for(int i=0; i<=getDim(); i++)
            if(vertices[i] != simpl.getConstVertices()[i]){
                return vertices[i] < simpl.getConstVertices()[i];
            }

    }

    return getDim() < simpl.getDim();
}

//--------------------------------------
//TopSimplex implementation

TopSimplex::TopSimplex()
{

}

TopSimplex::TopSimplex(vector<int> vertices)
{
    sort(vertices.begin(), vertices.end());
    this->vertices = vertices;
    adjacents = vector<int>(vertices.size(),INT_MAX-1);
}

TopSimplex::TopSimplex(int vertex)
{
    vertices.push_back(vertex);
}

int TopSimplex::getDimension()
{
    return vertices.size()-1;
}

int TopSimplex::get_nVertices()
{
    return vertices.size();
}

vector<int>& TopSimplex::getVertices()
{
    return vertices;
}


int TopSimplex::getAdjacent(int faceIndex) const
{
    return adjacents[faceIndex];
}

int TopSimplex::getVertexIndex(int vertex) const
{
    return vertices[vertex];
}

void TopSimplex::setAdjacent(int faceIndex, int simpl)
{
    adjacents[faceIndex]=simpl;
}




void TopSimplex::print_debug(){

    for(uint i=0; i<vertices.size(); i++)
        cout << vertices[i] << " ";
    cout << endl;

//    for(uint i=0; i<adjacents.size(); i++){
//        if(adjacents[i] < 0){
//            cout << "One adjacent: " << -(adjacents[i]+1) << endl;
//        }
//        else if(adjacents[i] == INT_MAX-1){
//            cout << "No adjacent" << endl;
//        }
//        else{
//            cout << "Multiple adjacents" << endl;
//        }
//    }

    cout << endl;
}
