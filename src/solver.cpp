#include "../include/solver.h"

#if FP_TRAP
# include <fenv.h>
#endif

solver::solver(
        const std::vector<objects::object *> &scene,
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

void solver::save(const std::string &prefix) {
    for (auto p : scene)
        p->save(prefix, _step);
    for (auto t : tracers) {
        t->walk(prefix, _step, gas(), [this] (vec p) -> const state {
            for (auto o : scene) {
                if (o->has_point(p))
                    return o->state_at(p);
            }
            state wrong;
            wrong.rho[0] = 1;
            return wrong;
        });
    }
}

double solver::estimate_timestep(const double dtlimit) {
    double dtmin = dtlimit;
    for (auto p : scene) {
        double dt = p->get_max_dt();
        if (dt < dtmin)
            dtmin = dt;
    }
    return cou * dtmin;
}

void solver::compute_flux(dir::Direction dir, const double dt) {
    for (auto p : scene) {
        p->compute_inner_flux(dir, dt / p->h(dir));
        p->compute_outer_flux(dir);
    }
}

void solver::compute_slope(dir::Direction dir) {
    for (auto p : scene)
        p->compute_inner_slope(dir);
}

void solver::integrate_by(dir::Direction dir, const double t, const double dt) {
    for (auto p : scene)
        p->integrate_by(dir, t, dt);
}

void solver::integrate_rhs(const double t, const double dt) {
    for (auto p : scene)
        p->integrate_rhs(t, dt);
}

void solver::integrate(const double dtlimit) {
    dt = estimate_timestep(dtlimit);
#if 0
    compute_slope(dir::X);
    compute_slope(dir::Y);
    compute_slope(dir::Z);
    compute_flux(dir::X, dt);
    compute_flux(dir::Y, dt);
    compute_flux(dir::Z, dt);
    integrate_by(dir::X, t, dt);
    integrate_by(dir::Y, t, dt);
    integrate_by(dir::Z, t, dt);
#else
    dir::Direction dirs[3];
    int v = step() % 6;
    dirs[0] = (v == 0 || v == 1) ? dir::X : ((v == 2 || v == 3) ? dir::Y : dir::Z);
    dirs[1] = (v == 2 || v == 4) ? dir::X : ((v == 0 || v == 5) ? dir::Y : dir::Z);
    dirs[2] = (v == 3 || v == 5) ? dir::X : ((v == 1 || v == 4) ? dir::Y : dir::Z);
/*    int v = rand() & 1;
    dirs[0] = v ? dir::Z : dir::X;
    dirs[1] = v ? dir::Y : dir::Y;
    dirs[2] = v ? dir::X : dir::Z; */
    for (int i = 0; i < 3; i++) {
        compute_slope(dirs[i]);
        compute_flux(dirs[i], dt);
        integrate_by(dirs[i], t, dt);
    }
#endif
    integrate_rhs(t, dt);
    t += dt;
    _step++;
}
