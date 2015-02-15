#include "../include/pipe.h"

namespace objects {

template struct pipe<NC>;

template<int nc>
void pipe<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, dir::Direction, double h, const double, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;

    cell.rhou(dir) -= dt * (right.fmom(dir) - left.fmom(dir)) / h;
    cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void pipe<nc>::integrate_rhs(state<nc> &cell, const state<nc> &source, const double t, const double dt) {
    object<nc>::integrate_rhs(cell, source, t, dt);

    double Dg = 4 * surface / perimeter;
    double Re = 1 + cell.rhou.norm() * Dg / this->gas().viscosity(cell);
    double cf = 0.0032 + 0.221 / pow(Re, 0.237);

    cf *= friction_coeff;

    cell.rhou(dir) -= dt * cell.density() * cell.velocity().norm2() * cf / 8;
}

template<int nc>
void pipe<nc>::compute_outer_flux(dir::Direction dir) {
    if (this->dir != dir)
        return;

    object<nc>::compute_outer_flux(dir);
}

template<int nc>
void pipe<nc>::integrate_by(dir::Direction dir, const double t, const double dt) {
    if (dir != this->dir)
        return;

    object<nc>::integrate_by(dir, t, dt);
}

template<int nc>
double pipe<nc>::get_max_dt() const {
    double cmax = 0;
    double vmax = 0;

    double hdir = this->h(dir);

    for (int i = 0; i < this->n(dir); i++) {
        const double c = this->gas().sound_speed(this->val(dir, i));
        const double v = fabs(this->val(dir, i).velocity()(dir));
        if (c > cmax)
            cmax = c;
        if (v > vmax)
            vmax = v;
    }

    return hdir / (cmax + vmax);
}

}
