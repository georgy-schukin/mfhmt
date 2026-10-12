#pragma once

#include "defs.h"
#include "async_op.h"

#include <list>
#include <tuple>
#include <memory>

class Reducer {
public:
    Reducer() {}

    void startNextReduction(const Double3 &data);
    Double3 finishTopReduction();
    void waitAll();

private:
    struct ReduceData {
        Double3 local;
        Double3 global;
        AsyncOp op;

        ReduceData(const Double3 &dt) : local(dt) {}
    };

private:
    std::list<std::shared_ptr<ReduceData>> reductions;
};
