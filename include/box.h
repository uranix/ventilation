#ifndef __BOX_H__
#define __BOX_H__

#include "vec.h"

#include <string>
#include <cassert>
#include <vector>

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

    ContactList *_side[3][2];

    vec ll, ur;
    vec h;

    bool closed[3][2];

    const std::string id;

    box(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
        : nx(nx), ny(ny), nz(nz), ll(ll), ur(ur), id(id)
    {
        h = (ur - ll) / vec(nx, ny, nz);

        assert(h.x > 0);
        assert(h.y > 0);
        assert(h.z > 0);

        for (int d = 0; d < 3; d++)
            for (int s = 0; s < 2; s++)
                closed[d][s] = false;

        for (int s = 0; s < 2; s++) {
            _side[0][s] = new ContactList[ny * nz];
            _side[1][s] = new ContactList[nx * nz];
            _side[2][s] = new ContactList[nx * ny];
        }
    }

    vec center(int i, int j, int k) {
        return ll + vec(i + .5, j + .5, k + .5) * h;
    }

    ContactList &side(int d, int s, int i, int j, int k) {
        if (d == 0) {
            return _side[d][s][j + ny * k];
        }
        if (d == 1) {
            return _side[d][s][i + nx * k];
        }
        return _side[d][s][i + nx * j];
    }

    ~box() {
        for (int d = 0; d < 3; d++)
            for (int s = 0; s < 2; s++)
                delete[] _side[d][s];
    }

    double get_min_h() const {
        return std::min(h.x, std::min(h.y, h.z));
    }
};

void connect(box &a, box &b);

#endif