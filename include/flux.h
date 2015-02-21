#ifndef __FLUX_H__
#define __FLUX_H__

#include "../include/state.h"
#include "../include/slope.h"

template<int nc>
struct predictor_flux {
    double fden[nc];
    vec fmom;
    double fener;

    void solve(const state<nc> &left, const state<nc> &right, const vec &norm, const gasinfo<nc> &gas);
};

template<int nc>
struct corrector_flux {
    double fden[nc];
    vec fmom;
    double fener;
    const double dt_h;
    corrector_flux(const double dt_h) : dt_h(dt_h) { }

    void solve(
            const state<nc> &left, const state<nc> &right,
            const slope<nc> &ls, const slope<nc> &cs, const slope<nc> &rs,
            const vec &norm, const gasinfo<nc> &gas);
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

    void add_outer(
            const state<nc> &left, const state<nc> &right,
            const vec &norm, const double Sfrac, const gasinfo<nc> &gas) 
    {
        outer_solve_add(left, right, norm, Sfrac, gas);
    }

    void set_inner(
            const state<nc> &left, const state<nc> &right,
            const slope<nc> &ls, const slope<nc> &cs, const slope<nc> &rs,
            const vec &norm, const double dt_h, const gasinfo<nc> &gas)
    {
        inner_solve_set(left, right, ls, cs, rs, norm, dt_h, gas);
    }

    void add_outer_reflect(
            optional<const state<nc> > _left, optional<const state<nc> > _right,
            const vec &norm, const double Sfrac, const gasinfo<nc> &gas)
    {
        assert((!_left) != (!_right));

        if (_left) {
            const state<nc> &left = *_left;
            state<nc> right(left);
            right.rhou = left.rhou - 2 * norm * dot(norm, left.rhou);
            outer_solve_add(left, right, norm, Sfrac, gas);
        } else {
            const state<nc> &right = *_right;
            state<nc> left(right);
            left.rhou = right.rhou - 2 * norm * dot(norm, right.rhou);
            outer_solve_add(left, right, norm, Sfrac, gas);
        }
    }

private:
    void outer_solve_add(const state<nc> &left, const state<nc> &right, const vec &norm, const double Sfrac, const gasinfo<nc> &gas) {
        predictor_flux<nc> pred;
        pred.solve(left, right, norm, gas);

        for (int i = 0; i < nc; i++)
            fdens[i] += Sfrac * pred.fden[i];

        fmom += Sfrac * pred.fmom;
        fener += Sfrac * pred.fener;
    }
    void inner_solve_set(
            const state<nc> &left, const state<nc> &right,
            const slope<nc> &ls, const slope<nc> &cs, const slope<nc> &rs,
            const vec &norm, const double dt_h, const gasinfo<nc> &gas)
    {
        predictor_flux<nc> pred;
        pred.solve(left, right, norm, gas);

        for (int i = 0; i < nc; i++)
            fdens[i] = pred.fden[i];

        fmom = pred.fmom;
        fener = pred.fener;

        corrector_flux<nc> corr(dt_h);
        corr.solve(left, right, ls, cs, rs, norm, gas);
        for (int i = 0; i < nc; i++)
            fdens[i] += corr.fden[i];

        fmom += corr.fmom;
        fener += corr.fener;
    }
};

#endif
