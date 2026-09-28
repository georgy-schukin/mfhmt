#pragma once

#include <vector>
#include <array>

class BlockDecomposition {
public:
    class Range {
    public:
        Range() {}
        Range(int s, int e):
            start(s), end(e) {}

        int size() const;
        int localStart(int g_index) const;
        int localEnd(int g_index) const;
        int toGlobal(int l_index) const;
        int toLocal(int g_index) const;
        bool hasIndex(int g_index) const;

    public:
        int start;
        int end;
    };

public:
    BlockDecomposition() {}
    BlockDecomposition(int size, int num_of_parts);
    BlockDecomposition(const std::vector<int> &sizes);

    int blockSize(int block_index) const;
    int blockShift(int block_index) const;
    int numOfBlocks() const;
    int fullSize() const;
    int localStart(int index, int block_index) const;
    int localEnd(int index, int block_index) const;
    int toGlobal(int index, int block_index) const;

    Range range(int block_index) const;

    BlockDecomposition grown(int amount) const;
    BlockDecomposition grown(int left, int right) const;
    BlockDecomposition multuplied(int amount) const;

private:
    static std::vector<int> computeSizes(int size, int num_of_parts);
    static std::vector<int> computeShifts(const std::vector<int> &sizes);

private:
    std::vector<int> sizes;
    std::vector<int> shifts;
};

class BlockDecomposition3D {
public:
    BlockDecomposition3D() {}
    BlockDecomposition3D(const BlockDecomposition &dx, const BlockDecomposition &dy, const BlockDecomposition &dz) :
        _decomps {dx, dy, dz} {
    }
    BlockDecomposition3D(int sx, int nx, int sy, int ny, int sz, int nz) :
        _decomps {BlockDecomposition(sx, nx), BlockDecomposition(sy, ny), BlockDecomposition(sz, nz)} {
    }

    const BlockDecomposition& decomp(int dim) const {
        return _decomps[dim];
    }

    int numOfBlocks(int dim) const {
        return decomp(dim).numOfBlocks();
    }

    int blockSize(int dim, int index) const {
        return decomp(dim).blockSize(index);
    }

    int blockShift(int dim, int index) const {
        return decomp(dim).blockShift(index);
    }

    BlockDecomposition::Range range(int dim, int index) const {
        return decomp(dim).range(index);
    }

private:
    std::array<BlockDecomposition, 3> _decomps;
};
