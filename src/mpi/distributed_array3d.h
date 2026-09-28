#pragma once

#include "defs.h"
#include "block_decomp.h"

#include <mpi.h>

#include <array>

class DistributedArray3D {
public:
    DistributedArray3D() {}
    DistributedArray3D(const BlockDecomposition3D &decomp, MPI_Comm cart_comm);

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

    const BlockDecomposition& decomp(int dim) const {
        return _decomp.decomp(dim);
    }

    const BlockDecomposition::Range& range(int dim) const {
        return _ranges[dim];
    }

    int rank() const {
        return _rank;
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

    void syncShadows(int dim);
    void syncShadows();

    DArray3 gather(int dst_rank) const;

private:
    int getShadowSize(const BlockDecomposition &dec) const {
        return dec.numOfBlocks() > 1 ? 1 : 0;
    }

private:
    DArray3 _data;
    BlockDecomposition3D _decomp;
    std::array<BlockDecomposition::Range, 3> _ranges;
    MPI_Comm _comm;
    int _rank;
    Index3 _index;
};
