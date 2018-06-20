#include "formangradient.h"

#include <phat/compute_persistence_pairs.h>

#include <phat/representations/vector_vector.h>
#include <phat/representations/vector_heap.h>
#include <phat/representations/vector_set.h>
#include <phat/representations/vector_list.h>
#include <phat/representations/sparse_pivot_column.h>
#include <phat/representations/heap_pivot_column.h>
#include <phat/representations/full_pivot_column.h>
#include <phat/representations/bit_tree_pivot_column.h>

#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>

#include <phat/helpers/dualize.h>

void FormanGradient::computePersistentHomology(){


    phat::boundary_matrix< phat::bit_tree_pivot_column > boundary_matrix;
    auto foo = bind(&FormanGradient::cmpSimplexesFiltr, this,_1,_2);

    list<implicitS> unsorted_set_cells;
    for(auto lvl : criticalS ){
        for(auto simpl : lvl.second){
            unsorted_set_cells.push_back(simpl);
        }
    }

    vector<implicitS> cells(unsorted_set_cells.begin(),unsorted_set_cells.end());
    sort(cells.begin(),cells.end(),foo);
    boundary_matrix.set_num_cols(cells.size());

    SUMap inverse_map = SUMap(foo);
    for(uint i=0; i<cells.size(); i++){
        inverse_map[cells[i]]=i;
    }

    //extract boundary matrices
    for(int i=0; i<cells.size(); i++){
        implicitS cell = cells[i];

        boundary_matrix.set_dim(i,cell.getDim());

        vector< phat::index > temp_col;
        if(cell.getDim() != 0){

            SSet boundary_cells;
            computeBoundaryCell(cell, boundary_cells);

            for(implicitS c : boundary_cells){
                assert(inverse_map.find(c) != inverse_map.end());
                temp_col.push_back(inverse_map[c]);
            }
        }

        sort(temp_col.begin(),temp_col.end());
        boundary_matrix.set_col(i,temp_col);
    }

    //run algorithm
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::twist_reduction >( pairs, boundary_matrix );
    pairs.sort();

    FILE* file = fopen("persistence_pairs.txt","w");
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        fprintf(file, "%f %f\n", simplexScalarValue(cells[pairs.get_pair( idx ).first]), simplexScalarValue(cells[pairs.get_pair( idx ).second]));
    }
    fclose(file);
}



