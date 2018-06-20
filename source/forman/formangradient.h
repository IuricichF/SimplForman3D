#ifndef FORMANGRADIENT_H
#define FORMANGRADIENT_H

#include <stack>
#include <math.h>

#include "gradientencoding.h"


typedef set<implicitS, boost::function<bool(const implicitS &, const implicitS &)>> SSet;
typedef map<implicitS, unsigned, boost::function<bool(const implicitS &, const implicitS &)>> SUMap;
typedef map<implicitS, bool, boost::function<bool(const implicitS &, const implicitS &)>> SBMap;


class FormanGradient
{
private:
    vector<uint> filtration; //for each vertex its filtration value

    GradientEncoding gradient;
    map<uint, SSet > criticalS;

    SimplicialComplex sc;


public:
    FormanGradient(int,char**);
    ~FormanGradient();

    //compute the Forman gradient on the input dataset, if batchTop=true the relation among vertices and incident top simplices
    //is computed once and stored explicitly (higher memory consumption but lower timings)
    void computeFormanGradient(bool batchTop);

    void computePersistentHomology();

    //print critical points and 1-skeleton
    void saveScalarFieldVTK(char* filename);
    void print1skeleton();

private:
    SSet* vertexLowerStar(uint vert, uint d); //compute the lower star of vert (only simplices of dimension=d)
    void getLowerStar(int v,SSet& lwStars); //split the lower star of a vertex v according to the sublevelsets of the function
    void homotopyExpansion(SSet&); //apply homotopy expansion on a set of simplices
    int numPairableLowerStar(const implicitS &next, const SSet& sset, implicitS& pair); //return the number of simplices pairable with next in the lower star, pair is one of these
    bool isPaired(const implicitS &simpl); //true if simpl is paired with another simplex
    void setPair(const implicitS &next, const implicitS &pair); //set the new gradient pair between next and pair (NOTE: next has to be bigger than pair)
    bool getPair(const implicitS &simpl, implicitS& next); //next is the simplex paired with simpl
    void freePair(const implicitS &next, const implicitS &pair); //remove pair (next,pair) from the gradient (NOTE: next has to be bigger than pair)

    uint simplexFiltration(const implicitS &simpl); //return the vector-valued filtration for a simplex. Each component is obtained as the maximum of the filtrations of its vertices
    float simplexScalarValue(const implicitS &simpl); //return the vector-valued function for a simplex. Each component is obtained as the maximum of the function values of its vertices

    void computeBoundaryCell(implicitS const& cell, SSet& descCells); //given a critical i-simplex return the critical (i-1)-cells connected to it

    //Output functions
    void computeDescendingCell(bool output,implicitS const& cell, SSet& desCells);
    void computeAscendingCell(bool output,implicitS const& cell, SSet& desCells);
    void print_out(string fileName, list<SSet> const& cells, int param, int dim);
    void out3cells(const list<SSet> &cells);
    void out2cells(const list<SSet> &cells, bool desc);
    void out1cells(const list<SSet> &cells, bool desc);
    void outCriticalPoints(const SSet &cells);
    void accurate_asc1cells(const list<SSet> &cells);

    //Simplex comparer - to organize simplices in maps and sets based on their function/filtration values
    bool cmpSimplexesFiltr(const implicitS& lhs, const implicitS& rhs); //true if the indexing of lhs is less than rhs (or if lhs is smaller than rhs). Used for homotopy expansion when working on a single sublevelset
    bool sortVerticesFiltration(const int& v1,const int& v2); //true if the filtration of v1 is less than v2
    bool filtrComparer(const pair<float,uint>& v1, const pair<float,uint>& v2) const; //true if v1 has a smaller value scalar value (float) or a smaller index (uint)


public:
    inline void outputInfos(){

        cout << endl << "-------Object infos-----------------------------------" << endl;
        cout << "The simplicial complex has: " << endl;
        cout << sc.getVerticesNum() << " vertices" << endl;
        for(auto index : sc.getTopSimplexesSet()){
            cout << sc.getTopSimplexesNum(index) << " " << index << "-simplices (only top simplices reported)" << endl;
        }
        cout << endl;

        cout << "The Forman gradient has critical cells: " << endl;
        for(auto lvl : criticalS){
            cout << "Dim: " << lvl.first << "  #: " << lvl.second.size() << endl;
        }
        cout << endl;
    }

    inline void outputCriticalPoints(){

        auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
        SSet points = SSet(foo);
        for(auto lvl : criticalS){
            points.insert(lvl.second.begin(),lvl.second.end());
        }
        outCriticalPoints(points);
    }

    inline void outputDescendingMorse(){

        for(auto lvl : criticalS){
            list<SSet> cells;
            for(auto s : lvl.second){
                if(lvl.first != 0){
                    SSet desc;
                    computeDescendingCell(true,s,desc);
                    cells.push_back(desc);
                }
            }

            switch (lvl.first) {
            case 1:
                out1cells(cells,true);
                break;
            case 2:
                out2cells(cells,true);
                break;
            case 3:
                out2cells(cells,true);
                break;
            default:
                break;
            }
        }
    }



};

#endif // FORMANGRADIENT_H
