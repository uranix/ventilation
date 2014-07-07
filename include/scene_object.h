#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "box.h"
#include "state.h"
#include "sloped_state.h"
#include "flux.h"

#include <fstream>

template<int nc>
struct functor {
    virtual void operator()(const vec &, state<nc> &) const = 0;
    virtual ~functor() { }
};

template<int nc>
class solver;

namespace objects {

template<int nc>
struct scene_object : public box {
    state<nc> *_oldstates;
    state<nc> *_states;
    state<nc> *_oldslopes[3];
    state<nc> *_slopes[3];
    state<nc> *_sources;
    flux<nc> *_fluxes[3];
    const solver<nc> *solver;

    scene_object(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
        : box(nx, ny, nz, ll, ur, id), solver(nullptr)
    {
        _states = new state<nc>[nx * ny * nz];
        _oldstates = new state<nc>[nx * ny * nz];
        for (int i = 0; i < 3; i++) {
            _slopes[i] = new state<nc>[nx * ny * nz];
            _oldslopes[i] = new state<nc>[nx * ny * nz];
        }
        _sources = new state<nc>[nx * ny * nz];

        _fluxes[0] = new flux<nc>[(nx + 1) * ny * nz];
        _fluxes[1] = new flux<nc>[(ny + 1) * nx * nz];
        _fluxes[2] = new flux<nc>[(nz + 1) * ny * nx];
    }

    void set_solver(const ::solver<nc> *solver) {
        this->solver = solver;
    }

    const vec &g() const {
        return solver->g();
    }

    const gasinfo<nc> &gas() const {
        return solver->gas();
    }

    virtual ~scene_object() {
        delete[] _states;
        delete[] _oldstates;
        delete[] _sources;

        for (int i = 0; i < 3; i++) {
            delete[] _slopes[i];
            delete[] _oldslopes[i];
            delete[] _fluxes[i];
        }
    }

    void fill(const functor<nc> &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++) {
                    sloped_state<nc> st = ref(i, j, k);
                    f(center(i, j, k), st.avg);
                }
    }

    void copy_explicit() {
        for (int i = 0; i < nx * ny * nz; i++) {
            _oldstates[i] = _states[i];
            _oldslopes[0][i] = _slopes[0][i];
            _oldslopes[1][i] = _slopes[1][i];
            _oldslopes[2][i] = _slopes[2][i];
        }
    }

    void average(const state<nc> &u0, state<nc> &u2) {
        for (int i = 0; i < nc; i++)
            u2.rho[i] = .5 * (u0.rho[i] + u2.rho[i]);
        u2.rhou = .5 * (u0.rhou + u2.rhou);
        u2.rhoE = .5 * (u0.rhoE + u2.rhoE);
    }

    void average_with_explicit() {
        for (int i = 0; i < nx * ny * nz; i++) {
            average(_oldstates[i], _states[i]);
            average(_oldslopes[0][i], _slopes[0][i]);
            average(_oldslopes[1][i], _slopes[1][i]);
            average(_oldslopes[2][i], _slopes[2][i]);
        }
    }

    void limit_slopes();

    void fill_sources(const functor<nc> &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), this->source(i, j, k));
    }

    const_sloped_state<nc> val(int i, int j, int k) const {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);

        int idx = i + (j + k * ny) * nx;
        return const_sloped_state<nc>(
                    _states[idx],
                    _slopes[0][idx],
                    _slopes[1][idx],
                    _slopes[2][idx]
                );
    }

    sloped_state<nc> ref(int i, int j, int k) {
        assert(i >= 0 && i < nx);
        assert(j >= 0 && j < ny);
        assert(k >= 0 && k < nz);

        int idx = i + (j + k * ny) * nx;
        return sloped_state<nc>(
                    _states[idx],
                    _slopes[0][idx],
                    _slopes[1][idx],
                    _slopes[2][idx]
                );
    }

    const_sloped_state<nc> val(dir::Direction dir, int i) const {
        if (dir == dir::X)
            return val(i, 0, 0);
        if (dir == dir::Y)
            return val(0, i, 0);
        return val(0, 0, i);
    }

    #define MAYBECONST
    #include "indexer.inc"
    #undef MAYBECONST
    #define MAYBECONST const
    #include "indexer.inc"
    #undef MAYBECONST

    virtual void compute_inner_fluxes();
    virtual void compute_outer_fluxes();

    virtual double get_max_dt() const;

    virtual void integrate(sloped_state<nc> cell, const flux<nc> &left, const flux<nc> &right,
            dir::Direction dir, double h, const double dt);
    virtual void integrate_rhs(sloped_state<nc> cell, const state<nc> &source, const double dt);
    virtual void integrate(const double dt);
    virtual void integrate_rhs(const double dt);

    template<typename T>
    void put(std::fstream &f, T value) const;

    void save(const std::string &prefix, const int step) const;

    void debug_avg() const;
};

#include "../src/scene_object.tpp"
#include "../src/vtk.tpp"

}

#endif
