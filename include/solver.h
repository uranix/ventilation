#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "object.h"
#include "tracer.h"

#include <vector>

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
    double time() const { return t; }
    double timestep() const { return dt; }
    int step() const { return _step; }

    double estimate_timestep(const double dtlimit);
    void integrate(const double dtlimit
            = objects::object<nc>::timestep_unconstrained);
    std::string version() const {
        #include "GitVersion.h"
        return VERSION;
    }
    void save(const std::string &prefix);
};

#endif
