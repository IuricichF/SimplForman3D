#ifndef TOPSIMPLEX_H
#define TOPSIMPLEX_H

#include <vector>
#include <iostream>
#include <map>
#include <climits>
#include <boost/dynamic_bitset.hpp>

using namespace std;

typedef unsigned int uint;

//explicit representation for a simplex, only simplexes that are explicitly encoded (vertices and top-simplexes) can be represented using this class
class explicitS{

private:
    int dim; //dimension of the top simplex
    int index; //index of the top simplex in the array of dim-simplexes

public:

    //constructors
    explicitS() : dim(-1), index(-1){}
    explicitS(int f, int s) : dim(f), index(s){}


    int getDim() const{ return dim;}
    //return: dimension of the simplex

    int getIndex() const{ return index;}
    //return: index of the simplex

    friend ostream& operator<<(ostream& out, const explicitS& simplex){
        out << "(" << simplex.getDim() << ", " << simplex.getIndex() << ") ";
        out.flush();
        return out;
    }
    //stream operator for printing out a simplex

    bool operator==(explicitS simpl){return dim == simpl.getDim() && index == simpl.getIndex();}
    //return: true if two explicitS represent the same top-simplex

    bool operator<(const explicitS simpl) const{
        if(dim == simpl.getDim())
            return index < simpl.getIndex();
        else{
            return dim < simpl.getDim();
        }
    }
    //operator used for orgnizing explicitS in sets and maps

    bool isValid() const{return dim!=-1;}
    //return true if the explicitS is valid
};


//implicit representation for a simplex based on a top-simplex. Every simplex in the simplicial complex can be represented with this class.

//Given a simplex s
//Given a top-simplex t, where either s is on the boundary of t or s == t
//the implicit represetation encodes the vertices of t composing s.

class implicitS{
private:
    vector<int> vertices; //vertices of the vertex

public:

    //constructors using:
    implicitS();
    implicitS(int v);
    implicitS(const vector<int>& simplexV);
    //top-simplex encoded as an explicitS and set of vertices of top forming the simplex

    uint getDim() const;
    //return: dimension of the simplex represented

    vector<int>* getVertices();
    //return: reference to the vertices

    vector<int> const& getConstVertices() const;
    //return: vertices encoded

    bool contains(const implicitS& simpl);

    bool operator==(const implicitS& simpl) const;
    //return: true if the two representations indicate the same simplex inside the top-simplex

    bool operator!=(const implicitS& simpl) const;
    //return: true if the two representations indicate the same simplex inside the top-simplex

    bool operator<(const implicitS& simpl) const;
//    //operator used to orgnanize implicitS in sets or maps

    friend ostream& operator<<(ostream& os, const implicitS& simplex)
    {
        vector<int> verts = simplex.getConstVertices();
        os << "[" << verts[0];
        for(int i=1; i<verts.size(); i++)
            os << "," << verts[i];
        os << "]";
        return os;
    }
    //stream operator for printing out a simplex
};


//top simplexes explicitly encoded in the data structure
class TopSimplex
{

private:
    vector<int> vertices; //set of vertex indices composing the top-simplex
    vector<int> adjacents; //set of adjacent top simplexes for each face, adjacents[i] contains an index to the top-simplexes adjacent on the ith face.

public:

    //constructors
    TopSimplex();
    TopSimplex(vector<int> vertices);
    TopSimplex(int vertex); // for vertices

    void setAdjacent(int fIndex, int top);
    //set a new adjacent top-simplex

    int getDimension();
    //return: dimension of the top-simplex

    int get_nVertices();
    //return: number of vertices in the top simplex

    vector<int>& getVertices();
    //return: set of vertex index composing the top-simplex

    int getAdjacent(int fIndex) const;
    //return: index of adjacent on face fIndex

    int getVertexIndex(int vertex) const;
    //return index of vertex inside the top-simplex

    void print_debug();
    //print out the top simplex

};

#endif // TOPSIMPLEX_H
