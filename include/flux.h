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
#if K_EPSILON_MODEL
    double frhok;
    double frhoeps;
    double gv2;
#endif

    flux() {
        zero();
    }

    void zero() {
        for (int i = 0; i < nc; i++)
            fdens[i] = 0;
        fmom = vec(0);
        fener = 0;
#if K_EPSILON_MODEL
        frhok = 0;
        frhoeps = 0;
        gv2 = 0;
#endif
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
            dir::Direction dir, const double dt_h, const double h, const gasinfo &gas)
    {
        inner_solve_set(left, right, ls, cs, rs, dir, dt_h, h, gas);
    }

    void add_outer_reflect(
            optional<const state > _left, optional<const state > _right,
            dir::Direction dir, const double Sfrac, const double h, const gasinfo &gas)
    {
        const vec norm(dir);
        assert((!_left) != (!_right));

        if (_left) {
            const state &left = *_left;
            state right(left);
            right.rhou = -left.rhou;
            outer_solve_add(left, right, dir, Sfrac, gas);
            add_k_epsilon(left, right, h, nullptr);
        } else {
            const state &right = *_right;
            state left(right);
            left.rhou = -right.rhou;
            outer_solve_add(left, right, dir, Sfrac, gas);
            add_k_epsilon(left, right, h, nullptr);
        }
    }

private:
    void add_k_epsilon(const state &left, const state &right, const double h, const double *pfden) {
#if K_EPSILON_MODEL
        const double sk = 1;
        const double se = 1.3;

        const auto &gv = (right.velocity() - left.velocity()) / h;
        gv2 = gv.norm2();
        const double mut = 0.5 * (left.mut() + right.mut());

        fmom -= mut * gv;
        double fdenstot = 0;
        if (pfden) {
            for (int i = 0; i < nc; i++)
                fdenstot += pfden[i];
        }

        frhok = fdenstot * ((fdenstot > 0) ? left.k : right.k);
        frhoeps = fdenstot * ((fdenstot > 0) ? left.eps : right.eps);
        frhok -= mut / sk * (right.k - left.k) / h;
        frhoeps -= mut / se * (right.eps - left.eps) / h;
#else
        (void)h;
#endif
    }
    void outer_solve_add(const state &left, const state &right, dir::Direction dir, const double Sfrac, const gasinfo &gas) {
        predictor_flux pred;
        slope cs(left, right, dir, gas);
        pred.solve(left, right, cs, dir, gas);

        for (int i = 0; i < nc; i++)
            fdens[i] += Sfrac * pred.fden[i];

        fmom += Sfrac * pred.fmom;
        fener += Sfrac * pred.fener;
#if K_EPSILON_MODEL
        double fdenstot = 0;
        for (int i = 0; i < nc; i++)
            fdenstot += pred.fden[i];

        frhok += Sfrac * fdenstot * ((fdenstot > 0) ? left.k : right.k);
        frhoeps += Sfrac * fdenstot * ((fdenstot > 0) ? left.eps : right.eps);
#endif
    }
    void inner_solve_set(
            const state &left, const state &right,
            const slope &ls, const slope &cs, const slope &rs,
            dir::Direction dir, const double dt_h, const double h, const gasinfo &gas)
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

        add_k_epsilon(left, right, h, pred.fden);
    }
};

#endif
