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

MPI_Datatype makeVectorType(MPI_Datatype src_type, int num_of_blocks, int block_size, MPI_Aint stride, int extent) {
    MPI_Datatype type, vtype;
    MPI_Type_create_hvector(num_of_blocks, block_size, stride, src_type, &type);
    if (extent > 0) {
        MPI_Type_create_resized(type, 0, extent, &vtype);
    } else {
        MPI_Type_dup(type, &vtype);
    }
    MPI_Type_commit(&vtype);
    return vtype;
}

MPI_Datatype makeVectorType(int num_of_blocks, int block_size, int stride, int extent) {
    return makeVectorType(MPI_DOUBLE, num_of_blocks, block_size, stride * sizeof(double), extent * sizeof(double));
}

MPI_Datatype makeBlockType(MPI_Datatype src_type, int block_size, int extent) {
    MPI_Datatype type, btype;
    MPI_Type_contiguous(block_size, src_type, &type);
    if (extent > 0) {
        MPI_Type_create_resized(type, 0, extent, &btype);
    } else {
        MPI_Type_dup(type, &btype);
    }
    MPI_Type_commit(&btype);
    return btype;
}

MPI_Datatype makeBlockType(int block_size, int extent) {
    return makeBlockType(MPI_DOUBLE, block_size, extent * sizeof(double));
}

MPI_Datatype makeSliceXType(const DArray3 &array) {
    // Slice in YZ plane
    return makeVectorType(array.size(1), array.size(2), array.fullSize(2), array.fullSize(1) * array.fullSize(2));
}

MPI_Datatype makeSliceYType(const DArray3 &array) {
    // Slice in XZ plane
    auto y_row_type = makeBlockType(array.size(2), array.fullSize(2));
    auto y_slice_type = makeVectorType(y_row_type, array.size(0), 1, array.fullSize(1) * array.fullSize(2) * sizeof(double), array.fullSize(2) * sizeof(double));
    MPI_Type_free(&y_row_type);
    return y_slice_type;
}

MPI_Datatype makeSliceZType(const DArray3 &array) {
    // Slice in XY plane
    auto z_col_type = makeVectorType(array.size(1), 1, array.fullSize(2));
    auto z_slice_type = makeVectorType(z_col_type, array.size(0), 1, array.fullSize(1) * array.fullSize(2) * sizeof(double), sizeof(double));
    MPI_Type_free(&z_col_type);
    return z_slice_type;
}

MPI_Datatype makeDataBlockType(const DArray3 &array) {
    auto x_slice_type = makeSliceXType(array);
    auto data_block_type = makeVectorType(x_slice_type, array.size(0), 1, array.fullSize(1) * array.fullSize(2) * sizeof(double), 0);
    MPI_Type_free(&x_slice_type);
    return data_block_type;
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


