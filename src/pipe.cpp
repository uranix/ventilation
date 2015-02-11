#include "pipe.h"

namespace objects {

template struct pipe<NC>;

template<int nc>
void pipe<nc>::compute_outer_fluxes() {
    const double tol = 1e-4;

    int i, j, k, idir;
    i = j = k = idir = 0;
    vec n = dir::to_vec(dir);
    double Sfrac = 1;

    this->flux_dir(dir, idir).zero();
    for (auto &z : this->side(dir, 0, i, j, k)) {
        Sfrac -= z.Sfrac;
        this->flux_dir(dir, idir).add(
                static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                this->val(dir, idir),
                n, z.Sfrac, this->gas(), this->g() * this->h(dir));
    }
    if (Sfrac > tol)
        this->flux_dir(dir, idir).add_reflect(
                this->val(dir, idir),
                false, n, Sfrac, this->gas(), this->g() * this->h(dir));

    idir = this->n(dir);
    dir::select(dir, i, j, k) = idir;
    Sfrac = 1;
    this->flux_dir(dir, idir).zero();
    for (auto &z : this->side(dir, 1, i, j, k)) {
        Sfrac -= z.Sfrac;
        this->flux_dir(dir, idir).add(
                this->val(dir, idir - 1),
                static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                n, z.Sfrac, this->gas(), this->g() * this->h(dir));
    }
    if (Sfrac > tol)
        this->flux_dir(dir, idir).add_reflect(
                this->val(dir, idir - 1),
                true, n, Sfrac, this->gas(), this->g() * this->h(dir));

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

template<int nc>
void pipe<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, dir::Direction, double h, const double, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;

    cell.rhou(dir) -= dt * (right.fmom(dir) - left.fmom(dir)) / h;
    cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void pipe<nc>::integrate(const double t, const double dt) {
    for (int i = 0; i < this->nx; i++)
        for (int j = 0; j < this->ny; j++)
            for (int k = 0; k < this->nz; k++) {
                if (dir == dir::X) integrate(this->ref(i, j, k), this->x_flux(i, j, k), this->x_flux(i+1, j, k), dir, this->h.x, t, dt);
                if (dir == dir::Y) integrate(this->ref(i, j, k), this->y_flux(i, j, k), this->y_flux(i, j+1, k), dir, this->h.y, t, dt);
                if (dir == dir::Z) integrate(this->ref(i, j, k), this->z_flux(i, j, k), this->z_flux(i, j, k+1), dir, this->h.z, t, dt);
            }
}

template<int nc>
void pipe<nc>::integrate_rhs(state<nc> &cell, const state<nc> &source, const double, const double dt) {
    cell.rhoE += dt * this->g().dot(cell.rhou);
    cell.rhou(dir) += dt * cell.density() * this->g()(dir);

    for (int i = 0; i < nc; i++)
        cell.rho[i] += dt * source.rho[i];
    cell.rhou += dt * source.rhou;
    cell.rhoE += dt * source.rhoE;

    double Dg = 4 * surface / perimeter;
    double Re = 1 + cell.rhou.norm() * Dg / this->gas().viscosity(cell);
    double cf = 0.0032 + 0.221 / pow(Re, 0.237);

    cf *= friction_coeff;

    cell.rhou(dir) -= dt * cell.density() * cell.velocity().norm2() * cf / 8;
}

}
