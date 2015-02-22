#ifndef __FLUX_H__
#define __FLUX_H__

#include "../include/state.h"
#include "../include/slope.h"

struct predictor_flux {
    double fden[nc];
    vec fmom;
    double fener;

    void solve(const state &left, const state &right,
            const slope &cs,
            dir::Direction dir, const gasinfo &gas);
};

struct corrector_flux {
    double fden[nc];
    vec fmom;
    double fener;
    const double dt_h;
    corrector_flux(const double dt_h) : dt_h(dt_h) { }

    void solve(
            const state &left, const state &right,
            const slope &ls, const slope &cs, const slope &rs,
            dir::Direction dir, const gasinfo &gas);
};

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
            const state &left, const state &right,
            dir::Direction dir, const double Sfrac, const gasinfo &gas)
    {
        outer_solve_add(left, right, dir, Sfrac, gas);
    }

    void set_inner(
            const state &left, const state &right,
            const slope &ls, const slope &cs, const slope &rs,
            dir::Direction dir, const double dt_h, const gasinfo &gas)
    {
        inner_solve_set(left, right, ls, cs, rs, dir, dt_h, gas);
    }

    void add_outer_reflect(
            optional<const state > _left, optional<const state > _right,
            dir::Direction dir, const double Sfrac, const gasinfo &gas)
    {
        const vec norm(dir);
        assert((!_left) != (!_right));

        if (_left) {
            const state &left = *_left;
            state right(left);
            right.rhou = left.rhou - 2 * norm * dot(norm, left.rhou);
            outer_solve_add(left, right, dir, Sfrac, gas);
        } else {
            const state &right = *_right;
            state left(right);
            left.rhou = right.rhou - 2 * norm * dot(norm, right.rhou);
            outer_solve_add(left, right, dir, Sfrac, gas);
        }
    }

private:
    void outer_solve_add(const state &left, const state &right, dir::Direction dir, const double Sfrac, const gasinfo &gas) {
        predictor_flux pred;
        slope cs(left, right, dir, gas);
        pred.solve(left, right, cs, dir, gas);

        for (int i = 0; i < nc; i++)
            fdens[i] += Sfrac * pred.fden[i];

        fmom += Sfrac * pred.fmom;
        fener += Sfrac * pred.fener;
    }
    void inner_solve_set(
            const state &left, const state &right,
            const slope &ls, const slope &cs, const slope &rs,
            dir::Direction dir, const double dt_h, const gasinfo &gas)
    {
        predictor_flux pred;
        pred.solve(left, right, cs, dir, gas);

        for (int i = 0; i < nc; i++)
            fdens[i] = pred.fden[i];

        fmom = pred.fmom;
        fener = pred.fener;

        corrector_flux corr(dt_h);
        corr.solve(left, right, ls, cs, rs, dir, gas);
        for (int i = 0; i < nc; i++)
            fdens[i] += corr.fden[i];

        fmom += corr.fmom;
        fener += corr.fener;
    }
};

#endif
