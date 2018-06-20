#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <map>

#include "topsimplex.h"

using namespace std;

class Vertex
{

private:
    map<int, vector<int> > partial_coboundaryTop; //set of top-simplexes, incident in the vertex, grouped by dimension. Only one simplex per connected component.
    vector<float> coordinates; //set of coordinates

public:

    //constructors
    Vertex();
    Vertex(vector<float>);

    void changeCoordinate(float,int);
    //set a new coordinate value

    vector<float>& getCoordinates();
    //return: vertex coordinates

    void addCoordinate(int, float);
    //return: add a new coordinate

    float getCoordinate(int) const;
    //return: ith coordinate

    int getCoboundaryTopNum();
    //return: number of top-simplexes, incident in the vertex, stored

    int getCoboundaryTopNum(int d);
    //return: number of top-simplexes of dimension d, incident in the vertex, stored

    int getCoboundaryMaxDim();
    //return: maximal simplex among the top-simplexes incident in the vertex

    void addPartialCoboundaryTop(int dim, int index);
    //insert a new top-simplex, of dimension dim, incident in the vertex

    map<int, vector<int> >& getPartialCoboundaryTopRelations();
    //return: top-simplexes in the vertex coboundary

    vector<int>* getPartialCoboundaryTop(int d);
    //return: top-simplexes of dimension d in the vertex coboundary

    float euclideanDistance(Vertex& v);
    //return: euclidean distance between two vertices

    void middlePoint(Vertex& v);

    bool operator<(const Vertex& v) const;

    friend ostream& operator<<(ostream& out, const Vertex& simplex){
        out << "(";
        for(auto f : simplex.coordinates)
            out << f << " ";
        out << ")";
        out.flush();
        return out;
    }
};

#endif // VERTEX_H
