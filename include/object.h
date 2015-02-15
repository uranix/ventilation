#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "box.h"
#include "state.h"
#include "flux.h"

#include <vector>
#include <algorithm>

template<int nc>
struct functor {
    virtual void operator()(const vec &, state<nc> &) const = 0;
    virtual ~functor() { }
};

template<int nc>
class solver;

namespace objects {

template<int nc>
struct object : public box {
private:
    std::vector<state<nc>> _states;
    std::vector<state<nc>> _sources;
    std::vector<flux<nc>> _fluxes[dir::DIR_END];
    const solver<nc> *slvr;
public:
    object(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id);
    virtual ~object() { }

    void set_solver(const ::solver<nc> *slvr);
    static constexpr double timestep_unconstrained = 1e20;
    const vec &g() const;
    const gasinfo<nc> &gas() const;

    void fill(const functor<nc> &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), ref(i, j, k));
    }

    void fill_sources(const functor<nc> &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), this->source(i, j, k));
    }

    const state<nc> &val(int i, int j, int k) const {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);

        int idx = i + (j + k * ny) * nx;
        return _states[idx];
    }

    state<nc> &ref(int i, int j, int k) {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);

        int idx = i + (j + k * ny) * nx;
        return _states[idx];
    }

    const state<nc> val(dir::Direction dir, int i) const {
        if (dir == dir::X)
            return val(i, 0, 0);
        if (dir == dir::Y)
            return val(0, i, 0);
        return val(0, 0, i);
    }

    const state<nc> &state_at(vec p) const {
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
    virtual void integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right,
            dir::Direction dir, double h, const double t, const double dt);
    virtual void integrate_rhs(state<nc> &cell, const state<nc> &source, const double t, const double dt);

    /* Per direction virtuals */
    virtual void compute_outer_flux(dir::Direction);
    virtual void integrate_by(dir::Direction dir, const double t, const double dt);

    /* Per object vitruals */
    virtual double get_max_dt() const;

    /* Regular methods */
    void compute_inner_flux(dir::Direction);
    void compute_inner_fluxes();
    void compute_outer_fluxes();

    void integrate(const double t, const double dt);
    void integrate_rhs(const double t, const double dt);

    void save(const std::string &prefix, const int step) const;
};

}

#endif
