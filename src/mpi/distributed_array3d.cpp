#include "distributed_array3d.h"
#include "comm_util.h"

DistributedArray3D::DistributedArray3D(const BlockDecomposition3D &decomp, const CartTopology<3> &tp) :
    _decomp(decomp),
    _topology(tp) {
    _index = _topology.thisCoord();
    const auto size_x = decomp.blockSize(0, _index[0]);
    const auto size_y = decomp.blockSize(1, _index[1]);
    const auto size_z = decomp.blockSize(2, _index[2]);
    const auto shadow_x = getShadowSize(decomp.decomp(0));
    const auto shadow_y = getShadowSize(decomp.decomp(1));
    const auto shadow_z = getShadowSize(decomp.decomp(2));
    _data = std::move(DArray3(size_x, size_y, size_z, shadow_x, shadow_y, shadow_z));
    _ranges = {decomp.range(0, _index[0]), decomp.range(1, _index[1]), decomp.range(2, _index[2])};
    _slice_types = {makeSliceXType(localArray()), makeSliceZType(localArray()), makeSliceZType(localArray())};
}

AsyncOps DistributedArray3D::syncShadows() {
    AsyncOps ops;
    ops.add(syncShadows(0));
    ops.add(syncShadows(1));
    ops.add(syncShadows(2));
    return ops;
}

AsyncOps DistributedArray3D::syncShadowsPrev() {
    AsyncOps ops;
    ops.add(syncShadowsPrev(0));
    ops.add(syncShadowsPrev(1));
    ops.add(syncShadowsPrev(2));
    return ops;
}

AsyncOps DistributedArray3D::syncShadowsNext() {
    AsyncOps ops;
    ops.add(syncShadowsNext(0));
    ops.add(syncShadowsNext(1));
    ops.add(syncShadowsNext(2));
    return ops;
}

AsyncOps DistributedArray3D::syncShadowsPrev(size_t dim) {
    AsyncOps ops;
    // Recv shadow from prev rank
    if (_topology.hasPrevNeighbor(dim)) {
        Index3 shadow_index {0, 0, 0};
        shadow_index[dim] = -1;
        auto op = recvSlice(_data, shadow_index, _slice_types[dim], _topology.prevNeighborRank(dim), _topology.mpiComm());
        ops.add(op);
    }
    // Send data to next rank
    if (_topology.hasNextNeighbor(dim)) {
        Index3 data_index {0, 0, 0};
        data_index[dim] = _data.size(dim) - 1;
        auto op = sendSlice(_data, data_index, _slice_types[dim], _topology.nextNeighborRank(dim), _topology.mpiComm());
        ops.add(op);
    }
    return ops;
}

AsyncOps DistributedArray3D::syncShadowsNext(size_t dim) {
    AsyncOps ops;
    // Recv shadow from next rank
    if (_topology.hasNextNeighbor(dim)) {
        Index3 shadow_index {0, 0, 0};
        shadow_index[dim] = _data.size(dim);
        auto op = recvSlice(_data, shadow_index, _slice_types[dim], _topology.nextNeighborRank(dim), _topology.mpiComm());
        ops.add(op);
    }
    // Send data to prev rank
    if (_topology.hasPrevNeighbor(dim)) {
        Index3 data_index {0, 0, 0};
        auto op = sendSlice(_data, data_index, _slice_types[dim], _topology.prevNeighborRank(dim), _topology.mpiComm());
        ops.add(op);
    }
    return ops;
}

AsyncOps DistributedArray3D::syncShadows(size_t dim) {
    AsyncOps ops;
    ops.add(syncShadowsPrev(dim));
    ops.add(syncShadowsNext(dim));
    return ops;
}

DArray3 DistributedArray3D::gather(int dst_rank) const {

}
