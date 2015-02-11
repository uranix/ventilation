#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "object.h"
#include "tracer.h"

#include <vector>
#if FP_TRAP
# include <fenv.h>
#endif

template<int nc>
class solver {
    std::vector<objects::object<nc> *> scene;
    std::vector<tracer *> tracers;
    const double cou;
    double t;
    double dt;
    int _step;

    gasinfo<nc> _gas;
    vec _g;
public:
    solver(const std::vector<objects::object<nc> *> &scene,
            const std::vector<tracer *> &tracers, const double cou);
    void set_gas(const gasinfo<nc> &gas) { this->_gas = gas; }
    void set_gravity(const vec &g) { this->_g = g; }
    const gasinfo<nc> &gas() const { return _gas; }
    const vec &g() const { return _g; }

    void compute_fluxes() {
        for (auto p : scene) {
            p->compute_inner_fluxes();
            p->compute_outer_fluxes();
        }
    }
    double estimate_timestep() {
        double dtmin = 1e20;
        for (auto p : scene) {
            double dt = p->get_max_dt();
            if (dt < dtmin)
                dtmin = dt;
        }
        return cou * dtmin;
    }
    void integrate() {
        compute_fluxes();
        dt = estimate_timestep();
        for (auto p : scene) {
            p->integrate_rhs(t, dt);
            p->integrate(t, dt);
        }
        t += dt;
        _step++;
    }
    std::string version() const {
        #include "GitVersion.h"
        return VERSION;
    }
    void save(const std::string &prefix);
    double time() const {
        return t;
    }
    double timestep() const {
        return dt;
    }
    int step() const {
        return _step;
    }
};

#endif
