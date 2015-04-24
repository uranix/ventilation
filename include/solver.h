#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "object.h"
#include "tracer.h"

#include <vector>

class solver {
    std::vector<objects::object *> scene;
    std::vector<tracer *> tracers;
    const double cou;
    double t;
    double dt;
    int _step;

    gasinfo _gas;
    vec _g;
public:
    solver(const std::vector<objects::object *> &scene,
            const std::vector<tracer *> &tracers, const double cou);
    void set_gas(const gasinfo &gas) { this->_gas = gas; }
    void set_gravity(const vec &g) { this->_g = g; }
    const gasinfo &gas() const { return _gas; }
    const vec &g() const { return _g; }
    double time() const { return t; }
    double timestep() const { return dt; }
    int step() const { return _step; }

    double estimate_timestep(const double dtlimit);
    void compute_slope(dir::Direction dir);
    void compute_flux(dir::Direction dir, const double dt);
    void special_flux(dir::Direction dir, const double dt);
    void integrate_by(dir::Direction dir, const double t, const double dt);
    void integrate_rhs(const double t, const double dt);

    void integrate(const double dtlimit
            = objects::object::timestep_unconstrained);
    std::string version() const {
        #include "GitVersion.h"
        return VERSION;
    }
    void save(const std::string &prefix);
};

#endif
