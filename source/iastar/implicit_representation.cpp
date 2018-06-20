#include "simplicialcomplex.h"

bool SimplicialComplex::areIncident(const implicitS &small, const implicitS &big){

    assert(small.getDim() < big.getDim());

    vector<int> result(big.getConstVertices().size());
    vector<int>::iterator it = set_intersection(small.getConstVertices().begin(), small.getConstVertices().end(),
                                                big.getConstVertices().begin(), big.getConstVertices().end(), result.begin());
    result.resize(it-result.begin());
    return result.size() == small.getConstVertices().size();
}

implicitS SimplicialComplex::toImplicit(const explicitS& simpl){
    if(simpl.getDim() == 0) return implicitS(simpl.getIndex());
    return implicitS( getTopSimplex(simpl).getVertices());
}



explicitS SimplicialComplex::toExplicit(const implicitS& simpl){

    vector<explicitS>* top = topStar(simpl);

    if(top->size()==1){
        return (*top)[0];
    }
    else{
        printf("The encoded simplex is neither a vertex nor a top simplex. Size of coboundary %d\n", top->size());
        return explicitS();
    }


}

bool SimplicialComplex::theSame(const implicitS& s1,const implicitS& s2){

    if(s1.getDim() == s2.getDim())
        return equal(s1.getConstVertices().begin(), s1.getConstVertices().end(),s2.getConstVertices().begin());
    return false;
}


vector<explicitS>* SimplicialComplex::topStar(const implicitS& simpl){

//    if(simpl.getDim() == 0){
//        return topStar(toExplicit(simpl));
//    }
//    else{
        vector<int> vertices = simpl.getConstVertices();

        if(topPerVertex.size() > 0){
            vector<explicitS> cobBig = topPerVertex[vertices[0]];

            for(int i=1; i<vertices.size(); i++){

                vector<explicitS> other = topPerVertex[vertices[i]];

                vector<explicitS> result(cobBig.size());
                vector<explicitS>::iterator it = set_intersection(cobBig.begin(), cobBig.end(),
                                 other.begin(), other.end(), result.begin());
                result.resize(it-result.begin());
                cobBig=result;
            }

            return new vector<explicitS>(cobBig.begin(), cobBig.end());
        }
        else{

            map<explicitS, uint> top_found;

            for(uint i=0; i<vertices.size(); i++){

                set<explicitS> top_per_vertex;

                for(int d=simpl.getDim()+1; d<=getVertex(vertices[i]).getCoboundaryMaxDim(); d++){

                    vector<explicitS>* top = topStar(explicitS(0,vertices[i]),d);
                    top_per_vertex.insert(top->begin(), top->end());

                    delete top;
                }

                for(set<explicitS>::iterator it=top_per_vertex.begin(); it!=top_per_vertex.end(); it++){
                    map<explicitS, uint>::iterator it2 = top_found.find(*it);
                    if(it2 == top_found.end())
                        top_found[*it]=1;
                    else
                        top_found[*it]=top_found[*it]+1;
                }

            }

            forward_list<explicitS> ret = forward_list<explicitS>();
            for(map<explicitS, uint>::iterator it = top_found.begin(); it!=top_found.end(); it++){
                if(it->second == vertices.size()){
                    ret.push_front(it->first);
                }
            }

            return new vector<explicitS>(ret.begin(), ret.end());
        }
    //}
}



vector<implicitS>* SimplicialComplex::coboundaryk(const implicitS& simplex,uint dim){


    vector<explicitS>* star = topStar(simplex);
    set<implicitS> simplexSet;
    for(int i=0; i<star->size(); i++){

        implicitS inTop(toImplicit((*star)[i]));

        vector<implicitS>* sset = getSubsets(inTop.getConstVertices(),dim+1);
        simplexSet.insert(sset->begin(),sset->end());

        for(auto cbd : *sset)
            if(cbd.contains(simplex))
                simplexSet.insert(cbd);

        delete sset;
    }
    delete star;

    return new vector<implicitS>(simplexSet.begin(), simplexSet.end());
}



vector<implicitS>* SimplicialComplex::boundaryk(const implicitS& simplex,uint dim){

    if(simplex.getDim() <= dim){
        printf("No simplexes of dimension %d on the boundary of a %d-simplex", dim, simplex.getDim());
        return new vector<implicitS>();
    }

    vector<implicitS>* sset = getSubsets(simplex.getConstVertices(),dim+1);

    return sset;
}


vector<implicitS>* SimplicialComplex::adjacents(const implicitS& simpl){

    set<implicitS> simplexesSet;

    if(simpl.getDim() > 0){
        vector<implicitS>* boundary = boundaryk(simpl,simpl.getDim()-1);

        for(vector<implicitS>::iterator it = boundary->begin(); it!=boundary->end(); it++){
            vector<implicitS>* cob = coboundaryk(*it,simpl.getDim());
            simplexesSet.insert(cob->begin(), cob->end());
            delete cob;
        }
        delete boundary;
    }
    else{
        vector<implicitS>* coboundary = coboundaryk(simpl,1);
        for(vector<implicitS>::iterator it = coboundary->begin(); it!=coboundary->end(); it++){
            vector<implicitS>* cob = boundaryk(*it,0);
            simplexesSet.insert(cob->begin(), cob->end());
            delete cob;
        }

        delete coboundary;
    }

    simplexesSet.erase(simpl);
    return new vector<implicitS>(simplexesSet.begin(), simplexesSet.end());
}



uint32_t SimplicialComplex::nextSubset(uint32_t v) {
    uint32_t t = v | (v - 1);
    return (t + 1) | (((~t & -~t) - 1) >> (__builtin_ctz(v) + 1));
}

vector<int> SimplicialComplex::filter(const vector<int>& v, uint32_t mask) {
    vector<int> res;
    while (mask) {
        res.push_back(v[__builtin_ctz(mask)]);
        mask &= mask - 1;
    }
    return res;
}

vector<implicitS>* SimplicialComplex::getSubsets(const vector<int>& arr, uint32_t k) {
    vector<implicitS>* s = new vector<implicitS>();
    uint32_t max = (1 << arr.size());
    for (uint32_t v = (1 << k) - 1; v < max; v = nextSubset(v)) {
        s->push_back(implicitS(filter(arr, v)));
    }
    return s;
}
