#pragma once

#include <mpi.h>

#include "defs.h"
#include "block_decomp.h"
#include "cart_topology.h"
#include "async_op.h"

MPI_Datatype makeVectorTypeT(MPI_Datatype src_type, int num_of_blocks, int block_size, MPI_Aint stride, int extent = -1);
MPI_Datatype makeVectorType(int num_of_blocks, int block_size, int stride, int extent = -1);
MPI_Datatype makeBlockTypeT(MPI_Datatype src_type, int block_size, int extent = -1);
MPI_Datatype makeBlockType(int block_size, int extent = -1);

MPI_Datatype makeSliceXType(const DArray3 &array);
MPI_Datatype makeSliceYType(const DArray3 &array);
MPI_Datatype makeSliceZType(const DArray3 &array);

MPI_Datatype makeDataBlockType(const DArray3 &array);
MPI_Datatype makeDataBlockType(int sx, int sy, int sz, int fx, int fy, int fz);

AsyncOp sendSlice(DArray3 &arr, const Index3 &src_index, MPI_Datatype slice_type, int neigh_rank, MPI_Comm comm);
AsyncOp recvSlice(DArray3 &arr, const Index3 &dst_index, MPI_Datatype slice_type, int neigh_rank, MPI_Comm comm);

DArray3 gatherArray(const DArray3 &local_data, const BlockDecomposition3D &decomp, const CartTopology<3> &tp, int root = 0);
