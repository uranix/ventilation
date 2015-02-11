#ifndef __BOX_H__
#define __BOX_H__

#include "vec.h"

#include <string>
#include <cassert>
#include <vector>
#include <iostream>

struct box {
    const int nx, ny, nz;

    struct Contact {
        struct box *other;
        float Sfrac;
        int ri, rj, rk;
        Contact() { }
        Contact(struct box *other, float Sfrac, int ri, int rj, int rk)
            : other(other), Sfrac(Sfrac), ri(ri), rj(rj), rk(rk)
        { }
    };

    int n(dir::Direction dir) const {
        if (dir == dir::X)
            return nx;
        if (dir == dir::Y)
            return ny;
        return nz;
    }

    typedef std::vector<Contact> ContactList;

    ContactList *_side[dir::DIR_END][dir::SIDE_END];

    vec ll, ur;
    vec h;

    bool closed[dir::DIR_END][dir::SIDE_END];

    const std::string id;

    box(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
        : nx(nx), ny(ny), nz(nz), ll(ll), ur(ur), id(id)
    {
        h = (ur - ll) / vec(nx, ny, nz);

        if (h.x <= 0 || h.y <= 0 || h.z <= 0) {
            std::cerr << "Object " << id << " has a negative dimension!" << std::endl;
            abort();
        }

        for (auto d : dir::DIRECTIONS)
            for (auto s : dir::SIDES)
                closed[d][s] = false;

        for (auto s : dir::SIDES) {
            _side[dir::X][s] = new ContactList[ny * nz];
            _side[dir::Y][s] = new ContactList[nx * nz];
            _side[dir::Z][s] = new ContactList[nx * ny];
        }
    }

    bool has_point(vec p) {
        if (p.x < ll.x || p.y < ll.y || p.z < ll.z)
            return false;
        if (p.x > ur.x || p.y > ur.y || p.z > ur.z)
            return false;
        return true;
    }

    void locate_point(vec p, int &i, int &j, int &k, vec &ofs) const {
        vec dp = p - ll;
        dp /= h;
        i = static_cast<int>(dp.x);
        j = static_cast<int>(dp.y);
        k = static_cast<int>(dp.z);
        if (i < 0) i = 0;
        if (j < 0) j = 0;
        if (k < 0) k = 0;
        if (i >= nx) i = nx - 1;
        if (j >= ny) j = ny - 1;
        if (k >= nz) k = nz - 1;
        dp -= vec(i + .5, j + .5, k + .5);
        dp *= 2;
        ofs = dp;
    }

    vec center(int i, int j, int k) {
        return ll + vec(i + .5, j + .5, k + .5) * h;
    }

    ContactList &side(dir::Direction d, dir::Side s, int i, int j, int k) {
        if (d == dir::X)
            return _side[d][s][j + ny * k];
        if (d == dir::Y)
            return _side[d][s][i + nx * k];
        return _side[d][s][i + nx * j];
    }

    ~box() {
        for (auto d : dir::DIRECTIONS)
            for (auto s : dir::SIDES)
                delete[] _side[d][s];
    }
};

void connect(box &a, box &b);

#endif
