#pragma once

#include "defs.h"
#include "block_decomp.h"
#include "cart_topology.h"
#include "async_op.h"

#include <mpi.h>

#include <array>

class DistributedArray3D {
public:
    //DistributedArray3D() {}
    DistributedArray3D(const BlockDecomposition3D &decomp, const CartTopology<3> &tp);

    double* data() {
        return _data.data();
    }

    const double* data() const {
        return _data.data();
    }

    size_t size(size_t dim) const {
        return _data.size(dim);
    }

    size_t shadowSize(size_t dim) const {
        return _data.shadowSize(dim);
    }

    size_t fullSize(size_t dim) const {
        return _data.fullSize(dim);
    }

    size_t size() const {
        return _data.size();
    }

    template <typename Index>
    size_t at(Index x, Index y, Index z) const {
        return _data.at<Index>(x, y, z);
    }

    double& operator[](size_t index) {
        return _data[index];
    }

    const double& operator[](size_t index) const {
        return _data[index];
    }

    double& operator()(const Index3 &ind) {
        return _data(ind[0], ind[1], ind[2]);
    }

    const double& operator()(const Index3 &ind) const {
        return _data(ind[0], ind[1], ind[2]);
    }

    template <typename Index>
    double& operator()(Index x, Index y, Index z) {
        return _data(x, y, z);
    }

    template <typename Index>
    const double& operator()(Index x, Index y, Index z) const {
        return _data(x, y, z);
    }

    typename std::vector<double>::iterator begin() {
        return _data.begin();
    }

    typename std::vector<double>::iterator end() {
        return _data.end();
    }

    const BlockDecomposition& decomp(size_t dim) const {
        return _decomp.decomp(dim);
    }

    const BlockDecomposition3D& decomp() const {
        return _decomp;
    }

    const BlockDecomposition::Range& range(size_t dim) const {
        return _ranges[dim];
    }

    int rank() const {
        return _topology.thisRank();
    }

    const Index3 &index() const {
        return _index;
    }

    const DArray3& localArray() const {
        return _data;
    }

    DArray3& localArray() {
        return _data;
    }

    bool hasIndex(int gx, int gy, int gz) const {
        return range(0).hasIndex(gx) && range(1).hasIndex(gy) && range(2).hasIndex(gz);
    }

    Index3 toLocal(int gx, int gy, int gz) const {
        return Index3 {range(0).toLocal(gx), range(1).toLocal(gy), range(2).toLocal(gz)};
    }

    AsyncOps syncShadows(size_t dim);
    AsyncOps syncShadowsPrev(size_t dim);
    AsyncOps syncShadowsNext(size_t dim);
    AsyncOps syncShadows();
    AsyncOps syncShadowsPrev();
    AsyncOps syncShadowsNext();

    DArray3 gather(int dst_rank) const;

private:
    int getShadowSize(const BlockDecomposition &dec) const {
        return dec.numOfBlocks() > 1 ? 1 : 0;
    }

private:
    DArray3 _data;
    const BlockDecomposition3D &_decomp;
    const CartTopology<3> &_topology;
    std::array<BlockDecomposition::Range, 3> _ranges;
    Index3 _index;
    std::array<MPI_Datatype, 3> _slice_types;
};
