#include "formangradient.h"


void FormanGradient::out3cells(const list<SSet> &cells){
    print_out("descending3cells.vtk",cells,10,4);
}

void FormanGradient::out2cells(const list<SSet> &cells, bool desc){
    if(desc)
        print_out("descending2cells.vtk",cells,5,3);
    else
        print_out("ascending2cells.vtk",cells,5,3);
}

void FormanGradient::out1cells(const list<SSet> &cells, bool desc){
    if(desc)
        print_out("descending1cells.vtk",cells,3,2);
    else{
        print_out("ascending1cells.vtk",cells,3,2);
        accurate_asc1cells(cells);
    }
}


void FormanGradient::outCriticalPoints(const SSet &cells){

    list<Vertex > coords;
    for(auto s : cells){
        vector<int> vertices = s.getConstVertices();
        Vertex v = sc.getVertex(vertices[0]);
        for(int i=0; i<vertices.size(); i++)
            v.middlePoint(sc.getVertex(vertices[i]));
        coords.push_back(v);
    }

    FILE* file = fopen("criticalpoints.vtk","w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", coords.size());

    for(auto c : coords){
        vector<float> coords = c.getCoordinates();
        fprintf(file, "%f %f %f \n", coords[0],coords[1],coords[2]);
    }

    fprintf(file, "\n\n");

    fprintf(file, "POINT_DATA %d \n", coords.size());
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "cp_type 1 %d float\n", coords.size());

    for(auto s : cells){
        fprintf(file, "%d ", s.getDim());
    }
    fprintf(file, "\n\n");


    fclose(file);
}


void FormanGradient::print_out(string fileName, list<SSet> const& cells, int param, int dim){

    map<int,int> oldToNew;
    vector<int> newToOld;

    int triangles = 0;

    int index=0;
    for(auto sets : cells){
        triangles+=sets.size();
        for(auto s : sets){
            vector<int> vertices = s.getConstVertices();
            for(auto v : vertices){

                if(oldToNew.find(v) == oldToNew.end())
                    oldToNew[v]=index++;
            }
        }
    }

    newToOld = vector<int>(oldToNew.size());
    for(auto v : oldToNew)
        newToOld[v.second]=v.first;


    int vertex_number=oldToNew.size();


    FILE* file = fopen(fileName.c_str(),"w");


    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertex_number);

    for(int i=0; i<vertex_number; i++){
        vector<float> coords = sc.getVertex(newToOld[i]).getCoordinates();
        fprintf(file, "%f %f %f \n", coords[0],coords[1],coords[2]);
    }

    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", triangles, triangles*(dim+1));

    for(auto sets : cells){
        for(auto c : sets){
            vector<int> vertices = c.getConstVertices();
            fprintf(file, "%d ", dim);
            for(int k=0; k<dim; k++){
               fprintf(file, "%d ", oldToNew[vertices[k]]);
            }
            fprintf(file, "\n");
        }
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", triangles);

    for(int i=0; i<triangles; i++)
        fprintf(file, "%d ", param);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", vertex_number);
    fprintf(file, "FIELD FieldData 1\n");

    fprintf(file, "originalField 1 %d float\n",vertex_number);

    for(int j=0; j<vertex_number; j++){
        fprintf(file, "%f ", sc.getVertex(newToOld[j]).getCoordinates().back());
    }



    fprintf(file, "\n\n");

    fprintf(file, "\n\n");

    fprintf(file, "CELL_DATA %d \n", triangles);
    fprintf(file, "FIELD FieldData 1\n");
    fprintf(file, "descending_2_cells 1 %d int \n", triangles);

    int reg=0;
    for(auto sets : cells){
        for(auto c : sets){
           fprintf(file, "%d ", reg);
        }
        reg++;
    }

    fclose(file);
}


void FormanGradient::accurate_asc1cells(const list<SSet> &cells){

    map<Vertex,int> vertices;
    list<pair<int,int> > edges;

    int index=0;
    for(auto sets : cells){

        for(auto s : sets){
            vector<explicitS>* tri_top = sc.topStar(s);
            Vertex v1,v2;
            if(tri_top->size() == 2){
                v1 = sc.barycenter((*tri_top)[0]);
                v2 = sc.barycenter((*tri_top)[1]);
            }
            else if(tri_top->size() == 1){
                v1 = sc.barycenter((*tri_top)[0]);
                v2 = sc.barycenter(s);
            }
            else{
                cout << "NO " << tri_top->size() << endl;
            }
            delete tri_top;

            int indexV1,indexV2;
            auto it1 = vertices.find(v1);
            if(it1 == vertices.end()){
                vertices[v1]=index;
                indexV1 = index++;
            }
            else
                indexV1 = it1->second;

            auto it2 = vertices.find(v2);
            if(it2 == vertices.end()){
                vertices[v2]=index;
                indexV2 = index++;
            }
            else
                indexV2 = it2->second;

            edges.push_back(pair<int,int>(indexV1,indexV2));
        }
    }

    vector<Vertex> vord(vertices.size());
    for(auto v : vertices)
        vord[v.second]=v.first;


    FILE* file = fopen("ascending1cells_correct.vtk","w");


    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vord.size());

    for(auto v : vord){
        fprintf(file, "%f %f %f \n", v.getCoordinate(0),v.getCoordinate(1),v.getCoordinate(2));
    }

    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", edges.size(), edges.size()*3);

    for(auto e : edges){
            fprintf(file, "2 %d %d\n", e.first, e.second);
    }
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", edges.size());

    for(int i=0; i<edges.size(); i++)
        fprintf(file, "3 ");
    fprintf(file, "\n\n");

    fclose(file);
}



void FormanGradient::saveScalarFieldVTK(char* filename){

    FILE* file = fopen(filename,"w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", sc.getVerticesNum());

    for(int i=0; i<sc.getVerticesNum(); i++){
        vector<float> coords = sc.getVertex(i).getCoordinates();
        fprintf(file, "%f %f %f \n", coords[0],coords[1],coords[2]);
    }

    fprintf(file, "\n\n");

    int sizeTop=0;
    int sizeIndices=0;
    for(auto d : sc.getTopSimplexesSet())
    {
        sizeTop += sc.getTopSimplexesNum(d);
        sizeIndices += sc.getTopSimplexesNum(d)*(d+1);
    }

    fprintf(file, "CELLS %d %d\n", sizeTop, sizeIndices);

    for(auto d : sc.getTopSimplexesSet())
    {
        for(auto s : sc.getTopSimplices(d)){
            vector<int> verts = s.getVertices();
            fprintf(file, "%d ", verts.size());
            for( int i=0; i<verts.size(); i++)
                fprintf(file, "%d ", verts[i]);
            fprintf(file,"\n");
        }
    }

    fprintf(file, "CELL_TYPES %d\n", sizeTop);
    for(auto d : sc.getTopSimplexesSet())
    {
        for(auto s : sc.getTopSimplices(d)){
            if(s.getDimension() == 1)
                fprintf(file, "3 ");
            else if(s.getDimension() == 2)
                fprintf(file, "5 ");
            else if(s.getDimension() == 3)
                fprintf(file, "10 ");
        }
    }
    fprintf(file, "\n");

    fprintf(file, "POINT_DATA %d \n", sc.getVerticesNum());
    fprintf(file, "FIELD FieldData 1\n");

    fprintf(file, "originalField 1 %d float\n",sc.getVerticesNum());

    for(int j=0; j<sc.getVerticesNum(); j++){
        fprintf(file, "%f ", sc.getVertex(j).getCoordinates().back());
    }
    fprintf(file,"\n");


    fclose(file);
}


void FormanGradient::print1skeleton(){

    list<SSet> asc;
    list<SSet> desc;

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    SSet criticals = SSet(foo);

    for(auto lvl :criticalS){
        for(auto s : lvl.second){

            criticals.insert(s);

            if(lvl.first == 1){
                SSet ascCells,descCells;
                computeAscendingCell(true,s,ascCells);
                computeDescendingCell(true,s,descCells);

                asc.push_back(ascCells);
                desc.push_back(descCells);
            }
        }
    }

    out1cells(desc, true);
    outCriticalPoints(criticals);
    accurate_asc1cells(asc);
}







