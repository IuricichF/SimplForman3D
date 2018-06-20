#include "simplicialcomplex.h"

SimplicialComplex::SimplicialComplex()
{

}

int SimplicialComplex::getComplexDim(){

    return getTopSimplexesSet().back();
}

int SimplicialComplex::getVerticesNum(){
    return vertices.size();
}

int SimplicialComplex::getTopSimplexesNum(int dim){
    if(realIndex.find(dim) == realIndex.end())
        return 0;
    else
        return topSimplexes[realIndex[dim]].size();
}

vector<int> SimplicialComplex::getTopSimplexesSet(){
    vector<int> result;
    for(auto r : realIndex){
        result.push_back(r.first);
    }
    return result;
}

void SimplicialComplex::buildDataStructure(){

    topPerVertex = vector<vector<explicitS> >();

    int adjRelStart=0;
    for(uint i=0; i<topSimplexes.size(); i++){

        int dim = topSimplexes[i][0].getDimension();

        //build adjacency relations
        list<forward_list<int> > adjacentRelationsLocal;
        vector<vector<int> > sorted_faces(topSimplexes[i].size()*(dim+1), vector<int>(dim+2));

        for(uint j=0; j<topSimplexes[i].size(); j++){

            TopSimplex tS = topSimplexes[i][j];
            vector<int> topVertices = tS.getVertices();
            for(int f=0; f<tS.getDimension()+1; f++){
                set<int> copy(topVertices.begin(), topVertices.end());
                copy.erase(topVertices[f]);
                vector<int> face(copy.begin(), copy.end());
                face.push_back(j);
                face.push_back(f);

                sorted_faces[(j*(dim+1))+f]=face;
            }
        }

        sort(sorted_faces.begin(), sorted_faces.end());

        uint j=1;
        while(j<sorted_faces.size()){

            vector<pair<int,int> > adjTopSimplexes;

            bool equal=true;
            do{
                adjTopSimplexes.push_back(pair<int,int>(sorted_faces[j-1][dim], sorted_faces[j-1][dim+1]));

                if(j == sorted_faces.size())
                    break;

                for(int k=0; k<dim; k++){
                    if(sorted_faces[j-1][k] != sorted_faces[j][k]){
                        equal=false;
                        break;
                    }
                }
                j++;


            } while(equal);


            if(adjTopSimplexes.size() == 2){
                getTopSimplex(explicitS(dim,adjTopSimplexes[0].first)).setAdjacent(adjTopSimplexes[0].second,-(adjTopSimplexes[1].first+1));
                getTopSimplex(explicitS(dim,adjTopSimplexes[1].first)).setAdjacent(adjTopSimplexes[1].second,-(adjTopSimplexes[0].first+1));
            }
            else if(adjTopSimplexes.size() > 2){
                forward_list<int> local_list;
                for(uint k=0; k<adjTopSimplexes.size(); k++){
                    local_list.push_front(adjTopSimplexes[k].first);
                    getTopSimplex(explicitS(dim,adjTopSimplexes[k].first)).setAdjacent(adjTopSimplexes[k].second,adjRelStart);
                }
                adjacentRelationsLocal.push_back(local_list);
                adjRelStart++;
            }
        }

        if(adjacentRelationsLocal.size() > 0)
            adjRelations.insert(adjRelations.end(), adjacentRelationsLocal.begin(), adjacentRelationsLocal.end());

        sorted_faces.clear();
        adjacentRelationsLocal.clear();

        vector<set<int> > incidentTop(vertices.size(),set<int>());
        //build partial incidence relations for vertices
        for(uint j=0; j<topSimplexes[i].size(); j++){

            TopSimplex tS = topSimplexes[i][j];
            for(int v=0; v<tS.getDimension()+1; v++){
                incidentTop[tS.getVertexIndex(v)].insert(j);
            }
        }


        for(uint j=0; j<vertices.size(); j++){
            while(incidentTop[j].size() >0){
                int first = *(incidentTop[j].begin());
                forward_list<explicitS>* local_cluster = incidentCluster(explicitS(0,j), explicitS(dim,first));
                incidentTop[j].erase(first);

                vertices[j].addPartialCoboundaryTop(dim,first);
                for(forward_list<explicitS>::iterator it=local_cluster->begin(); it!=local_cluster->end(); it++){
                    assert(it->getDim() == dim);
                    incidentTop[j].erase(it->getIndex());
                }
            }
        }
    }
}

void SimplicialComplex::build_top_simplex(set<uint> setR, set<uint> setP, set<uint> setX, const vector<set<uint>* >& arcs, map<int, list<TopSimplex>* >* top_simplexes_local){

    if(setP.size() == 0 && setX.size() == 0) {
        //set setR as maximal

        if(setR.size() > 1){
            vector<int> new_top(setR.begin(), setR.end());
            sort(new_top.begin(), new_top.end());
            TopSimplex top(new_top);

            if(top_simplexes_local->find(top.getDimension()) == top_simplexes_local->end()){
                (*top_simplexes_local)[top.getDimension()]=new list<TopSimplex>();
            }
            (*top_simplexes_local)[top.getDimension()]->push_back(top);
        }
    }
    else if(setP.size() > 0){
        //chose here the pivot vertex for P
        vector<uint> vectorP(setP.begin(), setP.end());

        uint pivot=vectorP[0];
        for(int i=0;i<vectorP.size();i++){
            if(arcs[vectorP[i]]->size() > arcs[pivot]->size())
                pivot=vectorP[i];
        }

        for(int i=0;i<vectorP.size();i++){
              uint ui=vectorP[i];
              if(arcs[ui]->find(pivot) == arcs[ui]->end()){

                setP.erase(ui);
                set<uint> adjui(arcs[ui]->begin(), arcs[ui]->end());

                set<uint> newR = setR;
                newR.insert(ui);

                set<uint> newP;
                set_intersection(setP.begin(),setP.end(),adjui.begin(),adjui.end(), std::inserter(newP,newP.begin()));

                set<uint> newX;
                set_intersection(setX.begin(),setX.end(),adjui.begin(),adjui.end(), std::inserter(newX,newX.begin()));

                build_top_simplex(newR, newP, newX, arcs, top_simplexes_local);

                setX.insert(ui);
            }
        }
    }
}

Vertex& SimplicialComplex::getVertex(int vertex){
    return vertices[vertex];
}

vector<TopSimplex>& SimplicialComplex::getTopSimplices(int dim){

    if(realIndex.find(dim) == realIndex.end())
        assert(false);
    else
        return topSimplexes[realIndex[dim]];

}

TopSimplex&  SimplicialComplex::getTopSimplex(explicitS simpl){
    assert(simpl.getDim() > 0);
    if(realIndex.find(simpl.getDim()) == realIndex.end())
        assert(false);
    return topSimplexes[realIndex[simpl.getDim()]][simpl.getIndex()];
}

vector<explicitS>* SimplicialComplex::topAdjacent(const explicitS& simpl, uint face_index){

    forward_list<explicitS> list;
    int adj = getTopSimplex(simpl).getAdjacent(face_index);

    if(adj == INT_MAX-1){
        return new vector<explicitS>();
    }
    else if(adj < 0){
        list.push_front(explicitS(simpl.getDim(),-(adj+1)));
    }
    else{
        for(forward_list<int>::iterator it=adjRelations[adj].begin(); it!=adjRelations[adj].end(); it++){
            list.push_front(explicitS(simpl.getDim(),*it));
        }
    }

    return new vector<explicitS>(list.begin(), list.end());
}

forward_list<explicitS>* SimplicialComplex::incidentCluster(explicitS vertex,explicitS topS){

    //check TOREMOVE
    //topsimplex must be incident to vertex

    forward_list<explicitS>* ret = new forward_list<explicitS>();
    set<explicitS> visited; //alternatively we can use a local vector. We avoid global properties for using parallelization

    queue<explicitS> adjacentSimplexes;
    adjacentSimplexes.push(topS);
    visited.insert(topS);
    while(!adjacentSimplexes.empty()){

        TopSimplex& top = getTopSimplex(adjacentSimplexes.front());


        for(int i=0; i<top.get_nVertices(); i++){
            if(top.getVertexIndex(i) != vertex.getIndex()){
                vector<explicitS>* adjs = topAdjacent(adjacentSimplexes.front(),i);
                assert(adjs != NULL);
                for(vector<explicitS>::iterator it=adjs->begin(); it!=adjs->end(); it++){
                    if(visited.find(*it) == visited.end()){
                        visited.insert(*it);
                        adjacentSimplexes.push(*it);
                    }
                }
                delete adjs;
            }
        }

        adjacentSimplexes.pop();

    }

    ret->insert_after(ret->before_begin(), visited.begin(), visited.end());
    return ret;
}



vector<explicitS>* SimplicialComplex::topStar(const explicitS& vertex){

    assert(vertex.getDim() == 0);


    if(topPerVertex.size() > vertex.getIndex() && topPerVertex[vertex.getIndex()].size() != 0){
        return new vector<explicitS>(topPerVertex[vertex.getIndex()]);
    }

    forward_list<explicitS> ret;
    Vertex v = getVertex(vertex.getIndex());

    map<int, vector<int> > partialCob = v.getPartialCoboundaryTopRelations();
    vector<explicitS> simplexes(v.getCoboundaryTopNum());

    int j=0;
    for(map<int, vector<int> >::iterator it = partialCob.begin(); it != partialCob.end(); it++){
            for(uint k=0; k<it->second.size(); k++){
                simplexes[j]=explicitS(it->first, it->second[k]);
                j++;
            }
    }

    for(uint i=0; i<simplexes.size(); i++){
        forward_list<explicitS>* adjs = incidentCluster(vertex,simplexes[i]);
        assert(adjs != NULL);
        {
            ret.insert_after(ret.before_begin(), adjs->begin(), adjs->end());
        }

        delete adjs;
    }

    return new vector<explicitS>(ret.begin(),ret.end());
}


vector<explicitS>* SimplicialComplex::topStar(const explicitS& vertex, int dimension){

    forward_list<explicitS> ret = forward_list<explicitS>();
    Vertex v = getVertex(vertex.getIndex());

    map<int, vector<int> > partialCob = v.getPartialCoboundaryTopRelations();
    vector<explicitS> simplexes(partialCob[dimension].size());

    for(uint k=0; k<partialCob[dimension].size(); k++)
        simplexes[k]=explicitS(dimension, (partialCob[dimension])[k]);

    for(uint i=0; i<simplexes.size(); i++){

        forward_list<explicitS>* adjs = incidentCluster(vertex,simplexes[i]);
        ret.insert_after(ret.before_begin(),adjs->begin(), adjs->end());
        delete adjs;
    }

    return new vector<explicitS>(ret.begin(), ret.end());
}


vector<implicitS>* SimplicialComplex::boundaryk(const explicitS& simplex,uint dim){

    if(simplex.getDim() < dim){
        printf("No simplexes of dimensions %d on a %d-simplex", dim, simplex.getDim());
        return new vector<implicitS>();
    }

    vector<implicitS>* sset = getSubsets(getTopSimplex(simplex).getVertices(),dim);

    return sset;
}



vector<implicitS>* SimplicialComplex::coboundaryk(const explicitS& simplex,uint dim){
    return coboundaryk(toImplicit(simplex), dim);
}


vector<implicitS>* SimplicialComplex::inTop(const explicitS &topSimpl, uint dim){

    vector<implicitS>* vec = new vector<implicitS>();
    if(dim == topSimpl.getDim()){
        vec->push_back(implicitS(getTopSimplex(topSimpl).getVertices()));
        return vec;
    }
    else if(dim > topSimpl.getDim())
        return new vector<implicitS>();
    else{

    }
}

void SimplicialComplex::storeFullStar(){

    topPerVertex=vector<vector<explicitS> >(getVerticesNum());
    for(int i=0; i<getVerticesNum(); i++){
        vector<explicitS>* tops = topStar(explicitS(0,i));
        sort(tops->begin(),tops->end());
        topPerVertex[i]=*tops;
        delete tops;
    }
}

void SimplicialComplex::emptyFullStar(){

    topPerVertex.clear();
}


