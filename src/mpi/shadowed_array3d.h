#pragma once

#include <vector>
#include <array>
#include <cstddef>

template <typename T>
class ShadowedArray3D {
public:
    ShadowedArray3D() {}
    ShadowedArray3D(size_t sx, size_t sy, size_t sz, size_t shadow_x = 0, size_t shadow_y = 0, size_t shadow_z = 0) :
        _size {sx, sy, sz},
        _shadow {shadow_x, shadow_y, shadow_z},
        _data(fullSize(), T {}) {
    }
    ShadowedArray3D(const std::array<size_t, 3> &sz, const std::array<size_t, 3> &shadow = {0, 0, 0}) :
        _size(sz),
        _shadow(shadow),
        _data(fullSize(), T {}) {
    }

    template <typename Arr>
    void copy(const Arr &src, int dx = 0, int dy = 0, int dz = 0) {
        for (int i = 0; i < src.size(0); i++)
            for (int j = 0; j < src.size(1); j++)
                for (int k = 0; k < src.size(2); k++)
                    _data[at(i + dx, j + dy, k + dz)] = src(i, j, k);
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
        return size(0) * size(1) * size(2);
    }

    size_t shadowSize(int dim) const {
        return _shadow[dim];
    }

    size_t fullSize(int dim) const {
        return _size[dim] + 2 * _shadow[dim];
    }

    size_t fullSize() const {
        return fullSize(0) * fullSize(1) * fullSize(2);
    }

    template <typename Index>
    size_t at(Index x, Index y, Index z) const {
        return (x + _shadow[0]) * _size[1] * _size[2] + (y + _shadow[1]) * _size[2] + z + _shadow[2];
    }

    template <typename Index>
    size_t atRaw(Index x, Index y, Index z) const {
        return (x) * _size[1] * _size[2] + (y) * _size[2] + z;
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

    template <typename Index>
    T& raw(Index x, Index y, Index z) {
        return _data[atRaw(x, y, z)];
    }

    template <typename Index>
    const T& raw(Index x, Index y, Index z) const {
        return _data[atRaw(x, y, z)];
    }

    typename std::vector<T>::iterator begin() {
        return _data.begin();
    }

    typename std::vector<T>::iterator end() {
        return _data.end();
    }

private:
    std::array<size_t, 3> _size;
    std::array<size_t, 3> _shadow;
    typename std::vector<T> _data;
};
