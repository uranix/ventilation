#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "scene_object.h"

#include <vector>

template<int nc>
class solver {
    std::vector<objects::scene_object<nc> *> scene;
    const double C;
    double t;
    int _step;

    gasinfo<nc> _gas;
    vec _g;
public:
    solver(const std::vector<objects::scene_object<nc> *> &scene, const double C) : scene(scene), C(C) 
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
        for (auto p : scene)
            p->compute_inner_fluxes();
        for (auto p : scene)
            p->compute_outer_fluxes();
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
    void integrate(const double dt) {
        for (auto p : scene)
            p->integrate(dt);
        for (auto p : scene)
            p->integrate_rhs(dt);

        t += dt;
        _step++;
    }
    void save(const std::string &prefix) {
        for (auto p : scene)
            p->save(prefix, _step);
    }
    double time() const {
        return t;
    }
    int step() const {
        return _step;
    }
};

#endif
