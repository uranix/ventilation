#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "box.h"
#include "state.h"
#include "gasinfo.h"
#include "slope.h"
#include "flux.h"

#include <vector>
#include <algorithm>

struct functor {
    virtual void operator()(const vec &, state &) const = 0;
    virtual ~functor() { }
};

class solver;

namespace objects {

struct object : public box {
private:
    std::vector<state> _states;
    std::vector<state> _sources;
    std::vector<flux> _fluxes[dir::DIR_END];
    std::vector<slope> _slopes[dir::DIR_END];
    const solver *slvr;
public:
    object(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id);
    virtual ~object() { }

    void set_solver(const ::solver *slvr);
    static constexpr double timestep_unconstrained = 1e20;
    const vec &g() const;
    const gasinfo &gas() const;

    void fill(const functor &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), ref(i, j, k));
    }

    void fill_sources(const functor &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), this->source(i, j, k));
    }

    const state &val(int i, int j, int k) const {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);

        int idx = i + (j + k * ny) * nx;
        return _states[idx];
    }

    state &ref(int i, int j, int k) {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);

        int idx = i + (j + k * ny) * nx;
        return _states[idx];
    }

    const state val(dir::Direction dir, int i) const {
        if (dir == dir::X)
            return val(i, 0, 0);
        if (dir == dir::Y)
            return val(0, i, 0);
        return val(0, 0, i);
    }

    const state &state_at(vec p) const {
        int i, j, k;
        vec ofs;
        locate_point(p, i, j, k, ofs);
        return val(i, j, k);
    }

    #define MAYBECONST
    #include "indexer.inc"
    #undef MAYBECONST
    #define MAYBECONST const
    #include "indexer.inc"
    #undef MAYBECONST

    /* Per cell/face virtuals */
    virtual void integrate(state &cell, const flux &left, const flux &right,
            dir::Direction dir, double h, const double t, const double dt);
    virtual void integrate_rhs(state &cell, const state &source, const double t, const double dt);

    /* Per direction virtuals */
    virtual void compute_outer_flux(dir::Direction);
    virtual void integrate_by(dir::Direction dir, const double t, const double dt);

    /* Per object vitruals */
    virtual double get_max_dt() const;

    /* Regular methods */
    void compute_inner_flux(dir::Direction, const double dt_h);
    void compute_inner_slope(dir::Direction);
    void compute_inner_slopes();
    void compute_inner_fluxes();
    void compute_outer_fluxes();

    void integrate(const double t, const double dt);
    void integrate_rhs(const double t, const double dt);

    void save(const std::string &prefix, const int step) const;
};

}

#endif
