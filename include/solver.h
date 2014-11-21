#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "scene_object.h"
#include "tracer.h"

#include <vector>

template<int nc>
class solver {
    std::vector<objects::scene_object<nc> *> scene;
    std::vector<tracer *> tracers;
    const double C;
    double t;
    double dt;
    int _step;

    gasinfo<nc> _gas;
    vec _g;
public:
    solver(const std::vector<objects::scene_object<nc> *> &scene,
            const std::vector<tracer *> &tracers, const double C)
        : scene(scene), tracers(tracers), C(C)
    {
        for (auto p : scene)
            p->set_solver(this);
        t = 0;
        _step = 0;
    }
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
        return C * dtmin;
    }
    void integrate() {
#if SECOND_ORDER
        for (auto p : scene)
            p->copy_explicit();
#endif

        compute_fluxes();
        dt = estimate_timestep();
        for (auto p : scene) {
            p->integrate_rhs(t, dt);
            p->integrate(t, dt);
            p->limit_slopes();
        }
        t += dt;

#if SECOND_ORDER
        compute_fluxes();
        for (auto p : scene) {
            p->integrate_rhs(t, dt);
            p->integrate(t, dt);
            p->limit_slopes();
        }

        for (auto p : scene)
            p->average_with_explicit();
#endif

        _step++;
    }
    std::string version() const {
        #include "GitVersion.h"
        return VERSION;
    }
    void save(const std::string &prefix) {
        for (auto p : scene)
            p->save(prefix, _step);
        for (auto t : tracers) {
            t->template walk<nc>(prefix, _step, gas(), [this] (vec p) {
                for (auto o : scene) {
                    if (o->has_point(p))
                        return o->state_at(p);
                }
                state<nc> wrong;
                wrong.rho[0] = 1;
                return wrong;
            });
        }
    }
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
