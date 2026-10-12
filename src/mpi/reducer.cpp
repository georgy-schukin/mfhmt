#include "reducer.h"

#include <mpi.h>

#include <stdexcept>

void Reducer::startNextReduction(const Double3 &data) {
    auto rd = std::make_shared<ReduceData>(data);
    MPI_Request req;
    MPI_Iallreduce(rd->local.data(), rd->global.data(), 3, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &req);
    rd->op = AsyncOp(req);
    reductions.push_back(rd);
}

Double3 Reducer::finishTopReduction() {
    if (reductions.empty()) {
        throw std::runtime_error("No reduction to finish");
    }
    auto &rd = reductions.front();
    rd->op.wait();
    const auto result = rd->global;
    reductions.pop_front();
    return result;
}

void Reducer::waitAll() {
    for (auto &rd: reductions) {
        rd->op.wait();
    }
    reductions.clear();
}
