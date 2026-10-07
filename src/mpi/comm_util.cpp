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

MPI_Datatype makeVectorTypeT(MPI_Datatype src_type, int num_of_blocks, int block_size, MPI_Aint stride, int extent) {
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
    return makeVectorTypeT(MPI_DOUBLE, num_of_blocks, block_size, stride * sizeof(double), extent * sizeof(double));
}

MPI_Datatype makeBlockTypeT(MPI_Datatype src_type, int block_size, int extent) {
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
    return makeBlockTypeT(MPI_DOUBLE, block_size, extent * sizeof(double));
}

MPI_Datatype makeSliceXType(const DArray3 &array) {
    // Slice in YZ plane
    return makeVectorType(array.size(1), array.size(2), array.fullSize(2), array.fullSize(1) * array.fullSize(2));
}

MPI_Datatype makeSliceYType(const DArray3 &array) {
    // Slice in XZ plane
    auto z_col_type = makeVectorType(array.size(2), 1, array.fullSize(1));
    auto y_slice_type = makeVectorTypeT(z_col_type, array.size(0), 1, array.fullSize(1) * array.fullSize(2) * sizeof(double));
    MPI_Type_free(&z_col_type);
    return y_slice_type;
}

MPI_Datatype makeSliceZType(const DArray3 &array) {
    // Slice in XY plane
    auto y_row_type = makeBlockType(array.size(2));
    auto z_slice_type = makeVectorTypeT(y_row_type, array.size(0), 1, array.fullSize(1) * array.fullSize(2) * sizeof(double));
    MPI_Type_free(&y_row_type);
    return z_slice_type;
}

MPI_Datatype makeDataBlockType(const DArray3 &array) {
    return makeDataBlockType(array.size(0), array.size(1), array.size(2), array.fullSize(0), array.fullSize(1), array.fullSize(2));
}

MPI_Datatype makeDataBlockType(int sx, int sy, int sz, int fx, int fy, int fz) {
    auto x_slice_type = makeVectorType(sy, sz, fz);
    auto data_block_type = makeVectorTypeT(x_slice_type, sx, 1, fy * fz * sizeof(double));
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

DArray3 gatherArray(const DArray3 &local_data, const BlockDecomposition3D &decomp, const CartTopology<3> &tp, int root) {
    AsyncOps ops;
    auto send_type = makeDataBlockType(local_data);
    MPI_Request sreq;
    MPI_Isend(&local_data(0, 0, 0), 1, send_type, root, TAG_GATHER, tp.mpiComm(), &sreq);
    MPI_Type_free(&send_type);
    ops.add(sreq);

    DArray3 data;
    if (tp.thisRank() == root) {
        const auto fx = decomp.decomp(0).fullSize();
        const auto fy = decomp.decomp(1).fullSize();
        const auto fz = decomp.decomp(2).fullSize();
        data = DArray3(fx, fy, fz);
        for (int i = 0; i < tp.dim(0); i++) {
            for (int j = 0; j < tp.dim(1); j++) {
                for (int k = 0; k < tp.dim(2); k++) {
                    MPI_Request rreq;
                    const auto src_sx = decomp.decomp(0).blockSize(i);
                    const auto src_sy = decomp.decomp(1).blockSize(j);
                    const auto src_sz = decomp.decomp(2).blockSize(k);
                    auto recv_type = makeDataBlockType(src_sx, src_sy, src_sz, fx, fy, fz);
                    const auto src_rank = tp.rankFromCoord({i, j, k});
                    const auto dst_x = decomp.decomp(0).blockShift(i);
                    const auto dst_y = decomp.decomp(1).blockShift(j);
                    const auto dst_z = decomp.decomp(2).blockShift(k);
                    MPI_Irecv(&data(dst_x, dst_y, dst_z), 1, recv_type, src_rank, TAG_GATHER, tp.mpiComm(), &rreq);
                    MPI_Type_free(&recv_type);
                    ops.add(rreq);
                }
            }
        }
    }
    ops.wait();
    return data;
}



