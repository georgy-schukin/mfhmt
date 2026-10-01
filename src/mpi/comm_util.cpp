#include "comm_util.h"

#include <exception>
#include <stdexcept>

namespace {

enum Tags: int {
    TAG_SLICE = 0,
    TAG_GATHER,
    TAG_SCATTER
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

AsyncOp sendSlice(DArray3 &arr, const Index3 &src_index, MPI_Datatype slice_type, int neigh_rank, MPI_Comm comm) {
    MPI_Request req;
    MPI_Isend(&arr(src_index[0], src_index[1], src_index[2]), 1, slice_type, neigh_rank, TAG_SLICE, comm, &req);
    return AsyncOp(req);
}

AsyncOp recvSlice(DArray3 &arr, const Index3 &dst_index, MPI_Datatype slice_type, int neigh_rank, MPI_Comm comm) {
    MPI_Request req;
    MPI_Irecv(&arr(dst_index[0], dst_index[1], dst_index[2]), 1, slice_type, neigh_rank, TAG_SLICE, comm, &req);
    return AsyncOp(req);
}

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


