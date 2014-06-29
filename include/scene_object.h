#ifndef __SCENE_OBJECT_H__
#define __SCENE_OBJECT_H__

#include "box.h"
#include "state.h"
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
    state<nc> *_states;
    state<nc> *_sources;
    flux<nc> *_fluxes[3];
    const solver<nc> *solver;

    scene_object(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id)
        : box(nx, ny, nz, ll, ur, id), solver(nullptr)
    {
        _states = new state<nc>[nx * ny * nz];
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
        delete[] _sources;

        for (int i = 0; i < 3; i++) {
            delete[] _fluxes[i];
        }
    }

    void fill(const functor<nc> &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), (*this)(i, j, k));
    }

    void fill_sources(const functor<nc> &f) {
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                for (int k = 0; k < nz; k++)
                    f(center(i, j, k), this->source(i, j, k));
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

    virtual void integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, const double dt);
    virtual void integrate_rhs(state<nc> &cell, const state<nc> &source, const double dt);
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
