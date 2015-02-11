#ifndef __FLUX_H__
#define __FLUX_H__

#include "state.h"

template<int nc>
struct solver_flux {
    double fden[nc];
    vec fmom;
    double fener;

    void solve(const state<nc> &left, const state<nc> &right, const vec &norm, const gasinfo<nc> &gas);
};

template<int nc>
struct flux {
    double fdens[nc];
    vec fmom;
    double fener;

    flux() {
        zero();
    }

    void zero() {
        for (int i = 0; i < nc; i++)
            fdens[i] = 0;
        fmom = vec(0);
        fener = 0;
    }

    void solve_and_add(const state<nc> &left, const state<nc> &right, const vec &norm, const double Sfrac, const gasinfo<nc> &gas) {
        solver_flux<nc> iface;
        iface.solve(left, right, norm, gas);

        for (int i = 0; i < nc; i++)
            fdens[i] += Sfrac * iface.fden[i];

        fmom += Sfrac * iface.fmom;
        fener += Sfrac * iface.fener;
    }

    void add(const state<nc> &left, const state<nc> &right, const vec &norm, const double Sfrac, const gasinfo<nc> &gas, const vec &) {
        solve_and_add(left, right, norm, Sfrac, gas);
    }

    void add_reflect(const state<nc> &inner, bool inner_is_left, const vec &norm, const double Sfrac, const gasinfo<nc> &gas, const vec &) {
        state<nc> left, right;

        if (inner_is_left) {
            left = inner;
            right = left;
            right.rhou = left.rhou - 2 * norm * dot(norm, left.rhou);
            solve_and_add(left, right, norm, Sfrac, gas);
        } else {
            right = inner;
            left = right;
            left.rhou = right.rhou - 2 * norm * dot(norm, right.rhou);
            solve_and_add(left, right, norm, Sfrac, gas);
        }
    }
};

#endif
