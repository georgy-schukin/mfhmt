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
}

void DistributedArray3D::syncShadows() {
    syncShadows(0);
    syncShadows(1);
    syncShadows(2);
}

void DistributedArray3D::syncShadows(size_t dim) {

}

DArray3 DistributedArray3D::gather(int dst_rank) const {

}
