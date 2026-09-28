#include "comm_util.h"

#include <exception>
#include <stdexcept>

namespace {

enum Tags: int {
    TAG_PREV = 0,
    TAG_NEXT = 1,
    TAG_GATHER = 2,
    TAG_SCATTER = 3
};

}

MPI_Datatype makeDataVectorType(MPI_Datatype src_type, int num_of_blocks, int block_size, int stride, int extent) {
    MPI_Datatype type, vtype;
    MPI_Type_vector(num_of_blocks, block_size, stride, src_type, &type);
    if (extent > 0) {
        MPI_Type_create_resized(type, 0, extent * sizeof(double), &vtype);
    } else {
        MPI_Type_dup(type, &vtype);
    }
    MPI_Type_commit(&vtype);
    return vtype;
}

MPI_Datatype makeDataVectorType(int num_of_blocks, int block_size, int stride, int extent) {
    return makeDataVectorType(MPI_DOUBLE, num_of_blocks, block_size, stride, extent);
}

MPI_Datatype makeSliceXType(const DArray3 &array) {

}

MPI_Datatype makeSliceYType(const DArray3 &array) {

}

MPI_Datatype makeSliceZType(const DArray3 &array) {

}

MPI_Datatype makeDataBlockType(const DArray3 &array) {

}

void syncShadowsXPrev(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank) {

}

void syncShadowsXNext(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank) {

}

void syncShadowsYPrev(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank) {

}

void syncShadowsYNext(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank) {

}

void syncShadowsZPrev(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank) {

}

void syncShadowsZNext(DArray3 &arr, MPI_Datatype slice_type, int rank, int neigh_rank) {

}

/*void syncShadowsColsPrev(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank > 0) {
        // Recv in shadow from prev.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), size_t(0)), 1, col_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &r1);
    }
    if (rank < size - 1) {
        // Send data column to next.
        MPI_Isend(&arr(size_t(0), arr.size(1) - 1), 1, col_type, rank + 1, TAG_PREV, MPI_COMM_WORLD, &r2);
    }
    if (rank > 0) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsColsNext(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank < size - 1) {
        // Recv in shadow from next.
        MPI_Irecv(&arr.raw(arr.shadowSize(0), arr.fullSize(1) - 1), 1, col_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &r1);
    }
    if (rank > 0) {
        // Send data column to prev.
        MPI_Isend(&arr(0, 0), 1, col_type, rank - 1, TAG_NEXT, MPI_COMM_WORLD, &r2);
    }
    if (rank < size - 1) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank > 0) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsRowsPrev(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank > 0) {
        // Recv in shadow from prev.
        MPI_Irecv(&arr.raw(size_t(0), arr.shadowSize(1)), 1, row_type, rank - 1, TAG_PREV, MPI_COMM_WORLD, &r1);
    }
    if (rank < size - 1) {
        // Send data row to next.
        MPI_Isend(&arr(arr.size(0) - 1, size_t(0)), 1, row_type, rank + 1, TAG_PREV, MPI_COMM_WORLD, &r2);
    }
    if (rank > 0) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsRowsNext(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    MPI_Request r1, r2;
    if (rank < size - 1) {
        // Recv in shadow from next.
        MPI_Irecv(&arr.raw(arr.fullSize(0) - 1, arr.shadowSize(1)), 1, row_type, rank + 1, TAG_NEXT, MPI_COMM_WORLD, &r1);
    }
    if (rank > 0) {
        // Send data row to prev.
        MPI_Isend(&arr(0, 0), 1, row_type, rank - 1, TAG_NEXT, MPI_COMM_WORLD, &r2);
    }
    if (rank < size - 1) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
    }
    if (rank > 0) {
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
    }
}

void syncShadowsCols(DArray2 &arr, MPI_Datatype col_type, int rank, int size) {
    syncShadowsColsPrev(arr, col_type, rank, size);
    syncShadowsColsNext(arr, col_type, rank, size);
}

void syncShadowsRows(DArray2 &arr, MPI_Datatype row_type, int rank, int size) {
    syncShadowsRowsPrev(arr, row_type, rank, size);
    syncShadowsRowsNext(arr, row_type, rank, size);
}*/

DArray3 gatherArray(const DArray3 &local_data, const BlockDecomposition3D &decomp, int rank, int size, int root) {

}

/*DArray2 gatherArrayCols(const DArray2 &local_data, const BlockDecomposition &cols_decomp, int rank, int size, int root) {
    std::vector<MPI_Request> reqs;
    if (rank == root) {
        reqs.resize(size + 1);
    } else {
        reqs.resize(1);
    }

    auto send_type = makeDataVectorType(local_data);
    MPI_Isend(&local_data(0, 0), 1, send_type, 0, TAG_GATHER, MPI_COMM_WORLD, &reqs[0]);
    MPI_Type_free(&send_type);

    DArray2 data;
    if (rank == root) {
        data = DArray2(local_data.size(0), cols_decomp.fullSize());
        for (int r = 0; r < size; r++) {
            auto recv_type = makeDataVectorType(data.size(0), cols_decomp.getBlockSize(r), data.fullSize(1));
            MPI_Irecv(&data(0, cols_decomp.getBlockShift(r)), 1, recv_type, r, TAG_GATHER, MPI_COMM_WORLD, &reqs[r + 1]);
            MPI_Type_free(&recv_type);
        }
    }
    MPI_Waitall(reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
    return data;
}*/


