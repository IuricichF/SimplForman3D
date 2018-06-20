#include "formangradient.h"

void FormanGradient::computeBoundaryCell(implicitS const& cell, SSet& descCells){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    descCells = SSet(foo);

    stack<implicitS> qu;
    qu.push(cell);

    while(!qu.empty()){

        implicitS top = qu.top();
        qu.pop();

        vector<implicitS>* bd = sc.boundaryk(top,top.getDim()-1);
        for(auto s : *bd){
            implicitS next;
            if(getPair(s,next)){
                assert(isPaired(next) && isPaired(s));
                if(next.getDim() == cell.getDim() && !sc.theSame(next,top)){
                    qu.push(next);
                }
            }
            else{
                if(!isPaired(s)){
                    if(descCells.find(s) == descCells.end())
                        descCells.insert(s);
                    else
                        descCells.erase(s);
                }
            }
        }
        delete bd;
    }
}

void FormanGradient::computeDescendingCell(bool output,implicitS const& cell, SSet& desCells){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    queue<implicitS> qu;

    if(output){
        desCells = SSet(foo);
        desCells.insert(cell);
    }

    qu.push(cell);
    while (!qu.empty()) {

        implicitS top = qu.front();
        qu.pop();

        vector<implicitS>* bd = sc.boundaryk(top,top.getDim()-1);

        for(auto s : *bd){
            implicitS next;
            if(getPair(s,next))
            {
                assert(isPaired(next) && isPaired(s));
                if(next.getDim() == cell.getDim() && !sc.theSame(next,top)){
                    qu.push(next);
                    if(output){
                        desCells.insert(next);
                    }
                }
            }
        }
        delete bd;
    }
}


void FormanGradient::computeAscendingCell(bool output,implicitS const& cell, SSet& ascCells){

    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);
    SSet cellsPath = SSet(foo);

    if(output){
        ascCells = SSet(foo);
        cellsPath.insert(cell);
    }

    queue<implicitS> qu;
    qu.push(cell);

    while(!qu.empty()){
        implicitS top = qu.front();
        qu.pop();

        vector<implicitS>* bd = sc.coboundaryk(top,top.getDim()+1);

        for(auto s : *bd){
            implicitS next;
            if(getPair(s,next)){
                if(next.getDim() == cell.getDim() && !sc.theSame(next,top)){
                    qu.push(next);
                    if(output){
                        cellsPath.insert(next);
                    }
                }
            }
            else{
                if(cell.getDim()==1 ){
                    ascCells.insert(cellsPath.begin(),cellsPath.end());
                    cellsPath.clear();
                }
            }
        }
        delete bd;
    }
}

