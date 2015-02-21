#include "../include/solver.h"

#if FP_TRAP
# include <fenv.h>
#endif

template class solver<NC>;

template<int nc>
solver<nc>::solver(
        const std::vector<objects::object<nc> *> &scene,
        const std::vector<tracer *> &tracers, const double cou)
: scene(scene), tracers(tracers), cou(cou)
{
#if FP_TRAP
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif
    for (auto p : scene)
        p->set_solver(this);
    t = 0;
    _step = 0;
}

template<int nc>
void solver<nc>::save(const std::string &prefix) {
    for (auto p : scene)
        p->save(prefix, _step);
    for (auto t : tracers) {
        t->template walk<nc>(prefix, _step, gas(), [this] (vec p) -> const state<nc> {
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

template<int nc>
double solver<nc>::estimate_timestep(const double dtlimit) {
    double dtmin = dtlimit;
    for (auto p : scene) {
        double dt = p->get_max_dt();
        if (dt < dtmin)
            dtmin = dt;
    }
    return cou * dtmin;
}

template<int nc>
void solver<nc>::compute_flux(dir::Direction dir, const double dt) {
    for (auto p : scene) {
        p->compute_inner_flux(dir, dt / p->h(dir));
        p->compute_outer_flux(dir);
    }
}

template<int nc>
void solver<nc>::compute_slope(dir::Direction dir) {
    for (auto p : scene)
        p->compute_inner_slope(dir);
}

template<int nc>
void solver<nc>::integrate_by(dir::Direction dir, const double t, const double dt) {
    for (auto p : scene)
        p->integrate_by(dir, t, dt);
}

template<int nc>
void solver<nc>::integrate_rhs(const double t, const double dt) {
    for (auto p : scene)
        p->integrate_rhs(t, dt);
}

template<int nc>
void solver<nc>::integrate(const double dtlimit) {
    dt = estimate_timestep(dtlimit);
#if 0
    compute_flux(dir::X);
    compute_flux(dir::Y);
    compute_flux(dir::Z);
    integrate_by(dir::X, t, dt);
    integrate_by(dir::Y, t, dt);
    integrate_by(dir::Z, t, dt);
#else
    compute_slope(dir::X);
    compute_flux(dir::X, dt);
    integrate_by(dir::X, t, dt);
    compute_slope(dir::Y);
    compute_flux(dir::Y, dt);
    integrate_by(dir::Y, t, dt);
    compute_slope(dir::Z);
    compute_flux(dir::Z, dt);
    integrate_by(dir::Z, t, dt);
#endif
    integrate_rhs(t, dt);
    t += dt;
    _step++;
}
