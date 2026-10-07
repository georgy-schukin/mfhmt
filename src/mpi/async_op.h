#pragma once

#include <memory>
#include <vector>
#include <iostream>

#include <mpi.h>

class AsyncOp {
public:
    class Request {
    public:
        Request(MPI_Request req) :
            request(req) {
        }

        ~Request() {
            wait();
        }

        void wait() {
            if (!finished) {                                
                MPI_Wait(&request, MPI_STATUS_IGNORE);                
                finished = true;
            }
        }

        void cancel() {
            if (!finished) {
                MPI_Cancel(&request);
                finished = true;
            }
        }
    private:
        MPI_Request request;
        bool finished = false;
    };

public:
    AsyncOp() {}
    AsyncOp(MPI_Request req) {
        req_handle = std::make_shared<Request>(req);
    }

    void wait() {
        if (req_handle) {
            req_handle->wait();
            req_handle = nullptr;
        }
    }

    void cancel() {
        if (req_handle) {
            req_handle->cancel();
            req_handle = nullptr;
        }
    }
private:
    std::shared_ptr<Request> req_handle;
};

class AsyncOps {
public:
    AsyncOps() {}

    void add(const AsyncOp &op) {
        _ops.push_back(op);
    }

    void add(const AsyncOps &ops) {
        if (ops.size() > 0) {
            _ops.insert(_ops.end(), ops._ops.begin(), ops._ops.end());
        }
    }

    void add(MPI_Request req) {
        add(AsyncOp(req));
    }

    void wait() {        
        for (auto &op: _ops) {
            op.wait();
        }
        _ops.clear();
    }

    void cancel() {
        for (auto &op: _ops) {
            op.cancel();
        }
        _ops.clear();
    }

    size_t size() const {
        return _ops.size();
    }

private:
    std::vector<AsyncOp> _ops;
};
