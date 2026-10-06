#pragma once

#include <cstddef>
#include <array>

#include <mpi.h>

template <size_t Dims>
class CartTopology {
public:
    CartTopology(int rank, int size, const std::array<int, Dims> &dims = {}) :
        _rank(rank), _dims(dims) {
        std::array<int, Dims> periods {};
        MPI_Dims_create(size, static_cast<int>(Dims), _dims.data());
        MPI_Cart_create(MPI_COMM_WORLD, static_cast<int>(Dims), _dims.data(), periods.data(), 0, &_cart_comm);
        MPI_Cart_coords(_cart_comm, _rank, static_cast<int>(Dims), _coord.data());
    }

    int numOfDims() const {
        return Dims;
    }

    int dim(size_t dm) const {
        return _dims[dm];
    }

    const std::array<int, Dims>& thisCoord() const {
        return _coord;
    }

    int thisRank() const {
        return _rank;
    }

    int numOfNodes() const {
        int size = 0;
        MPI_Comm_size(_cart_comm, &size);
        return size;
    }

    bool hasPrevNeighbor(size_t dm) const {
        return _coord[dm] > 0;
    }

    bool hasNextNeighbor(size_t dm) const {
        return _coord[dm] + 1 < _dims[dm];
    }

    int prevNeighborRank(size_t dm) const {
        auto neigh_coord = _coord;
        neigh_coord[dm]--;
        return rankFromCoord(neigh_coord);
    }

    int nextNeighborRank(size_t dm) const {
        auto neigh_coord = _coord;
        neigh_coord[dm]++;
        return rankFromCoord(neigh_coord);
    }

    int rankFromCoord(const std::array<int, Dims> &coord) const {
        int rank;
        MPI_Cart_rank(_cart_comm, coord.data(), &rank);
        return rank;
    }

    MPI_Comm mpiComm() const {
        return _cart_comm;
    }

private:
    int _rank;
    std::array<int, Dims> _dims;
    std::array<int, Dims> _coord;
    MPI_Comm _cart_comm;
};
