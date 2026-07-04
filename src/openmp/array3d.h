#pragma once

#include <vector>
#include <array>
#include <cstddef>

template <typename T>
class Array3D {
public:
    Array3D() {}
    Array3D(size_t sx, size_t sy, size_t sz, const T &value = T {}) :
        _size {sx, sy, sz},
        _data(sx * sy * sz, value) {
    }
    Array3D(const std::array<size_t, 3> &sz, const T &value = T {}) :
        _size(sz),
        _data(sz[0] * sz[1] * sz[2], value) {
    }

    void populate(const T* raw_data, size_t data_sz) {
        for (size_t i = 0; i < data_sz; i++) {
            _data[i] = raw_data[i];
        }
    }

    T* data() {
        return _data.data();
    }

    const T* data() const {
        return _data.data();
    }

    size_t size(size_t dim) const {
        return _size[dim];
    }

    size_t size() const {
        return _data.size();
    }

    template <typename Index>
    size_t at(Index x, Index y, Index z) const {
        return x * _size[1] * _size[2] + y * _size[2] + z;
    }

    T& operator[](size_t index) {
        return _data[index];
    }

    const T& operator[](size_t index) const {
        return _data[index];
    }

    template <typename Index>
    T& operator()(Index x, Index y, Index z) {
        return _data[at(x, y, z)];
    }

    template <typename Index>
    const T& operator()(Index x, Index y, Index z) const {
        return _data[at(x, y, z)];
    }

    typename std::vector<T>::iterator begin() {
        return _data.begin();
    }

    typename std::vector<T>::iterator end() {
        return _data.end();
    }

private:
    std::array<size_t, 3> _size;
    typename std::vector<T> _data;
};
