#include "simplicialcomplex.h"

void SimplicialComplex::readOFF(const char* file){
    //reader for the OFF file

    realIndex.clear();

    ifstream fStream (file);
      if (fStream.is_open())
      {
        string line;
        int nV, nT, x;

        getline(fStream,line);
        if(line.compare(0,3,"OFF") != 0){
            cout << line << endl;
            printf("Wrong file: this is not an OFF file\n");
            exit(0);
        }

        //read number of vertices and number of top simplexes
        fStream >> nV;
        fStream >> nT;
        fStream >> x;

        vertices = vector<Vertex>(nV);
        getline(fStream,line);
        //read coordinates for each vertex

        for(int i=0; i<nV; i++){

            getline(fStream,line);
            istringstream iss(line);
            string coord;
            list<float> coords;
            while( getline(iss, coord, ' ')){
                coords.push_back(atof(coord.c_str()));
            }

            vector<float> coordinates(coords.begin(),coords.end());
            vertices[i] = Vertex(coordinates);
        }

        map<int, list<TopSimplex>* > topSimplexeslists=map<int, list<TopSimplex>* >();

        //read top simplexes
        int nVIndexes;
        for(int i=0; i<nT; i++){
            fStream >> nVIndexes;

            vector<int> topVertices(nVIndexes);

            for(int k=0; k<nVIndexes; k++){
                fStream >> topVertices[k];
            }

            TopSimplex topS(topVertices);

            if(topSimplexeslists.find(topS.getDimension()) == topSimplexeslists.end()){
                topSimplexeslists[topS.getDimension()]=new list<TopSimplex>();
            }
            topSimplexeslists[topS.getDimension()]->push_back(topS);

        }

        int dim=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            dim++;
        }


        fStream.close();
      }


      buildDataStructure();
      cout << "Complex built" << endl;
}

void SimplicialComplex::readTS(const char* file){
    //reader for the OFF file


    ifstream fStream (file);
      if (fStream.is_open())
      {
        string line;
        int nV, nT, x;

        //read number of vertices and number of top simplexes
        fStream >> nV;
        fStream >> nT;

        vertices = vector<Vertex>(nV);
        getline(fStream,line);
        //read coordinates for each vertex
        for(int i=0; i<nV; i++){

            getline(fStream,line);
            istringstream iss(line);
            string coord;
            list<float> coords;
            while( getline(iss, coord, ' ')){
                coords.push_back(atof(coord.c_str()));
            }

            vertices[i] = Vertex(vector<float>(coords.begin(),coords.end()));

        }

        //read top simplexes
        realIndex[3]=0;
        topSimplexes.push_back(vector<TopSimplex>(nT));

        for(int i=0; i<nT; i++){

            vector<int> topVertices(4);

            for(int k=0; k<4; k++){
                fStream >> topVertices[k];
            }

            TopSimplex topS(topVertices);
            topSimplexes[0][i]=topS;

        }

        fStream.close();
      }

      cout << "Building structure " << endl;
      buildDataStructure();
}


void SimplicialComplex::readGMV(const char* file){
    //reader for the GMV file


    ifstream fStream (file);
      if (fStream.is_open())
      {
        string line;
        getline(fStream,line);
        getline(fStream,line);
        getline(fStream,line);
        getline(fStream,line);
        getline(fStream,line);

        int nV, nT;
        fStream >> line;
        fStream >> nV;

        vertices = vector<Vertex>(nV);
        float coord;
        //read x coordinates for each vertex
        for(int i=0; i<nV; i++){
            vector<float> coords(3);
            fStream >> coords[0];
            vertices[i] = Vertex(coords);
        }
        //read y coordinates for each vertex
        for(int i=0; i<nV; i++){
            fStream >> coord;
            getVertex(i).changeCoordinate(coord,1);
        }
        //read z coordinates for each vertex
        for(int i=0; i<nV; i++){
            fStream >> coord;
            getVertex(i).changeCoordinate(coord,2);
        }

        fStream >> line;
        fStream >> nT;

        map<int, list<TopSimplex>* > topSimplexeslists;

        //read top simplexes
        int nVIndexes;
        for(int i=0; i<nT; i++){
            fStream >> line;
            fStream >> nVIndexes;

            vector<int> topVertices(nVIndexes);

            for(int k=0; k<nVIndexes; k++){
                fStream >> topVertices[k];
                topVertices[k] = topVertices[k]-1;
            }

            TopSimplex topS(topVertices);
            if(topSimplexeslists.find(topS.getDimension()) == topSimplexeslists.end()){
                topSimplexeslists[topS.getDimension()]=new list<TopSimplex>();
            }
            topSimplexeslists[topS.getDimension()]->push_back(topS);
        }

        int dim=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            dim++;
        }

        fStream.close();
      }

      buildDataStructure();
}


void SimplicialComplex::readPoints(char* name_in_file, float threshold, uint dimension){

    fstream fStream(name_in_file, ios::in);

    if (fStream.is_open())
    {
        list<Vertex >* points=new list<Vertex >();

        string line;
        double x;
        while( getline ( fStream , line ) )
        {
          vector<float> point;
          std::istringstream iss( line );
          while(iss >> x) { point.push_back(x); }
          points->push_back(Vertex(point));
        }
        fStream.close();

        vertices = vector<Vertex>(points->begin(), points->end());

        delete points;


        vector<set<uint>* > arcs(vertices.size());
        set<uint> ps;

        for(uint i=0; i<vertices.size(); i++){
            ps.insert(i);
            arcs[i] = new set<uint>();
        }


        int arcs_count=0;


        for(uint i=0; i<vertices.size(); i++){
            for(uint j=i+1; j<vertices.size(); j++){
                if(threshold >= vertices[i].euclideanDistance(vertices[j])){

                    arcs[i]->insert(j);

                    arcs[j]->insert(i);

                    arcs_count++;
                }
            }
        }


        map<int, list<TopSimplex>* > topSimplexeslists;
        build_top_simplex(set<uint>(), ps, set<uint>(), arcs, &topSimplexeslists);


        int dim=0;
        int totTop=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            totTop+=topSimplexes[dim].size();
            dim++;
        }

        buildDataStructure();
    }
}

void SimplicialComplex::readGraph(char* name_in_file){

    fstream fStream(name_in_file, ios::in);

    if (fStream.is_open())
    {
        string line;
        int vertN;
        int edges;
        getline ( fStream , line );

        std::istringstream iss( line );
        iss >> vertN;
        iss >> edges;

        vertices = vector<Vertex>(vertN);

        vector<set<uint>* > arcs(vertices.size());

        set<uint> points;
        for(uint i=0; i<vertices.size(); i++){
            arcs[i] = new set<uint>();
            points.insert(i);
        }


        map<int, int> vertMap;
        int v1,v2;
        int vReal=0;
        for(int i=0; i<edges; i++){
            getline ( fStream , line );

            std::istringstream iss( line );
            iss >> v1;
            iss >> v2;

            if(vertMap.find(v1) == vertMap.end()){
                vertMap[v1]=vReal;
                v1=vReal;
                vReal++;
            }
            else
                v1 = vertMap[v1];

            if(vertMap.find(v2) == vertMap.end()){
                vertMap[v2]=vReal;
                v2=vReal;
                vReal++;
            }
            else
                v2 = vertMap[v2];


            arcs[v1]->insert(v2);
            arcs[v2]->insert(v1);
        }


        map<int, list<TopSimplex>* > topSimplexeslists;
        build_top_simplex(set<uint>(), points, set<uint>(), arcs, &topSimplexeslists);

        int dim=0;
        int totTop=0;
        topSimplexes= vector<vector<TopSimplex> >(topSimplexeslists.size(), vector<TopSimplex>());
        for(map<int, list<TopSimplex>*>::iterator it=topSimplexeslists.begin(); it!=topSimplexeslists.end(); it++){

            realIndex[it->first]=dim;
            topSimplexes[dim]=vector<TopSimplex>(it->second->begin(), it->second->end());
            delete it->second;
            totTop+=topSimplexes[dim].size();
            dim++;
        }

        buildDataStructure();
    }
}

void SimplicialComplex::readIA(char* file){
    //reader for the OFF file

    ifstream fStream (file);
    if (fStream.is_open())
    {
        //number of vertices and number of different types of top simplexes
        int nV, nT;
        fStream >> nV;
        fStream >> nT;

        //for each type of top simplex, its real index in the vector of top simplexes
        int index, real;
        for(int i=0; i<nT; i++){
            fStream >> index;
            fStream >> real;
            realIndex[index]=real;
        }


        //for each vertex
        vertices = vector<Vertex>(nV);
        int coordDim;
        int cobDim, cobSize, cobInt, cobIndex;
        for(int i=0; i<nV; i++){

            //number of coordinates and values
            fStream >> coordDim;
            vector<float> coords(coordDim);
            for(int j=0; j<coordDim; j++)
                fStream >> coords[j];

            vertices[i] = Vertex(coords);
            fStream >> cobDim;
            for(int j=0; j<cobDim; j++){
                fStream >> cobInt;
                fStream >> cobSize;
                for(int k=0; k<cobSize; k++){
                    fStream >> cobIndex;
                    vertices[i].addPartialCoboundaryTop(cobInt,cobIndex);
                }
            }
        }

        //for each topsimplex dimension
        int topSize;
        int topDim;
        for(int i=0; i<nT; i++){
            //topSimplexes = vector<vector<TopSimplex> >(nT,vector<TopSimplex>());
            fStream >> topSize;
            vector<TopSimplex> tops(topSize);
            for(int j=0; j<topSize; j++){
                fStream >> topDim;
                vector<int> verts(topDim);
                for(int k=0; k<topDim; k++){
                    fStream >> verts[k];
                }
                tops[j]=TopSimplex(verts);
                int adjs;
                for(int k=0; k<topDim; k++){
                    fStream >> adjs;
                    tops[j].setAdjacent(k,adjs);
                }
            }
            topSimplexes.push_back(tops);
        }

        //for each non-manifold adjacent relation
        int adjRels, dim, indexes;
        fStream >> adjRels;
        adjRelations = vector<forward_list<int> >(adjRels);
        for(int i=0; i<adjRels; i++){
            forward_list<int> adj;
            fStream >> dim;
            for(int j=0; j<dim; j++){
                fStream >> indexes;
                adj.push_front(indexes);
            }
            adjRelations[i]=adj;
        }

        fStream.close();
    }
}


void SimplicialComplex::readSquareGrid(const char* file, int xRes, int yRes){

    int res = yRes*xRes;
    FILE* fp = fopen(file, "r");
    char* data = (char*) malloc (res*sizeof(char));
    fread((char *) data, res*sizeof(char), 1, fp);

    int vIndex=0;
    vertices = vector<Vertex>(res);

    realIndex[2]=0;
    topSimplexes.push_back(vector<TopSimplex>((yRes-1)*(xRes-1)*2));

    for(int j=0; j<yRes; j++){
        for(int i=0; i<xRes; i++){
            vector<float> coords = {(float)i, (float)j, (float)data[i + j*xRes], (float)data[i + j*xRes]};
            vertices[vIndex] = Vertex(coords);
            vIndex++;
        }
    }

    int tIndex=0;
    for(int j=0; j<yRes-1; j++){
        for(int i=0; i<xRes-1; i++){

            int v0 = i + j*xRes;

            int v1= (i+1) + j*xRes;
            int v2= i + (j+1)*xRes;
            int v3= (i+1) + (j+1)*xRes;

            vector<int> t0={v3, v2, v1};
            vector<int> t1={v0, v2, v1};

            topSimplexes[0][tIndex++]=t0;
            topSimplexes[0][tIndex++]=t1;
        }
    }

    buildDataStructure();
}

void SimplicialComplex::readCubicalGrid(const char *file, int xRes, int yRes, int zRes){

    int res = zRes*yRes*xRes;
    FILE* fp = fopen(file, "r");
    char* data = (char*) malloc (res*sizeof(char));
    fread((char *) data, res*sizeof(char), 1, fp);
    cout << "File read" << endl;
    int vIndex=0;
    vertices = vector<Vertex>(res);

    realIndex[3]=0;
    topSimplexes.push_back(vector<TopSimplex>((zRes-1)*(yRes-1)*(xRes-1)*6));

    for(int k=0; k<zRes; k++){
        for(int j=0; j<yRes; j++){
            for(int i=0; i<xRes; i++){
                vector<float> coords = {(float)i, (float)j, (float)k, (float)data[i + j*xRes + k*xRes*yRes]};
                vertices[vIndex] = Vertex(coords);
                vIndex++;
            }
        }
    }

    int tIndex=0;
    for(int k=0; k<zRes-1; k++){
        for(int j=0; j<yRes-1; j++){
            for(int i=0; i<xRes-1; i++){
                int v0 = i + j*xRes + k*xRes*yRes;

                int v1= (i+1) + j*xRes + k*xRes*yRes;
                int v2= i + (j+1)*xRes + k*xRes*yRes;
                int v3= (i+1) + (j+1)*xRes + k*xRes*yRes;;
                int v4= i + j*xRes + (k+1)*xRes*yRes;
                int v5= (i+1) + j*xRes + (k+1)*xRes*yRes;
                int v6= i + (j+1)*xRes + (k+1)*xRes*yRes;
                int v7= (i+1) + (j+1)*xRes + (k+1)*xRes*yRes;

                vector<int> t0={v0, v2, v6, v7};
                vector<int> t1={v7, v4, v6, v0};
                vector<int> t2={v0, v4, v5, v7};
                vector<int> t3={v7, v1, v5, v0};
                vector<int> t4={v7, v2, v3, v0};
                vector<int> t5={v0, v1, v3, v7};

                topSimplexes[0][tIndex++]=t0;
                topSimplexes[0][tIndex++]=t1;
                topSimplexes[0][tIndex++]=t2;
                topSimplexes[0][tIndex++]=t3;
                topSimplexes[0][tIndex++]=t4;
                topSimplexes[0][tIndex++]=t5;
            }
        }
    }

    buildDataStructure();
}


void SimplicialComplex::saveVTK(char* filename){

    FILE* file = fopen(filename,"w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", getVerticesNum());

    for(int i=0; i<getVerticesNum(); i++){
        vector<float> coords = getVertex(i).getCoordinates();
        fprintf(file, "%f %f %f \n", coords[0],coords[1],coords[2]);
    }

    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", getTopSimplexesNum(2), getTopSimplexesNum(2)*4);

    for(auto t : getTopSimplices(2)){
        vector<int> verts = t.getVertices();
        fprintf(file, "3 %d %d %d\n", verts[0],verts[1],verts[2]);
    }

    fprintf(file, "CELL_TYPES %d\n", getTopSimplexesNum(2));

    for(int i=0; i<getTopSimplexesNum(2); i++)
        fprintf(file, "5 ");
    fprintf(file, "\n\n");


    int nFields = getVertex(0).getCoordinates().size()-3;

    fprintf(file, "POINT_DATA %d \n", getVerticesNum());
    fprintf(file, "FIELD FieldData %d\n", nFields);

    for(int i=0; i<nFields; i++){
        fprintf(file, "originalField_%d 1 %d float\n",i+1,getVerticesNum());

        for(int j=0; j<getVerticesNum(); j++){
            fprintf(file, "%f ", getVertex(j).getCoordinates()[i+3]);
        }
        fprintf(file,"\n");
    }


    fclose(file);
}



