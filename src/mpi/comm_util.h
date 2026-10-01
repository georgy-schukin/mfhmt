#pragma once

#include <mpi.h>

#include "defs.h"
#include "block_decomp.h"

MPI_Datatype makeDataVectorType(MPI_Datatype src_type, int num_of_blocks, int block_size, int stride, int extent = -1);
MPI_Datatype makeDataVectorType(int num_of_blocks, int block_size, int stride, int extent = -1);
MPI_Datatype makeSliceXType(const DArray3 &array);
MPI_Datatype makeSliceYType(const DArray3 &array);
MPI_Datatype makeSliceZType(const DArray3 &array);
MPI_Datatype makeDataBlockType(const DArray3 &array);

void syncShadowsXPrev(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank);
void syncShadowsXNext(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank);
void syncShadowsYPrev(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank);
void syncShadowsYNext(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank);
void syncShadowsZPrev(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank);
void syncShadowsZNext(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank);

DArray3 gatherArray(const DArray3 &local_data, const BlockDecomposition3D &decomp, int rank, int size, int root = 0);

