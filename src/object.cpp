#include "../include/object.h"
#include "../include/solver.h"

#include <iostream>

namespace objects {

template struct object<NC>;

template<int nc>
object<nc>:: object(
        int nx, int ny, int nz,
        const vec &ll, const vec &ur, const std::string &id)
:
    box(nx, ny, nz, ll, ur, id), slvr(nullptr)
{
    std::cout << "Using scene object `" << id << "' with dims = "
        << "(" << nx << ", " << ny << ", " << nz << ") and h = "
        << "(" << h.x << ", " << h.y << ", " << h.z << ")" << std::endl;

    _states.resize(nx * ny * nz);
    _sources.resize(nx * ny * nz);

    _fluxes[dir::X].resize((nx + 1) * ny * nz);
    _fluxes[dir::Y].resize((ny + 1) * nx * nz);
    _fluxes[dir::Z].resize((nz + 1) * ny * nx);
}

template<int nc>
void object<nc>::set_solver(const ::solver<nc> *slvr) { this->slvr = slvr; }

template<int nc>
const vec &object<nc>::g() const { return slvr->g(); }

template<int nc>
const gasinfo<nc> &object<nc>::gas() const { return slvr->gas(); }

/* Per cell/face virtuals */
template<int nc>
void object<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, dir::Direction, double h, const double, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
    cell.rhou -= dt * (right.fmom - left.fmom) / h;
    cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void object<nc>::integrate_rhs(state<nc> &cell, const state<nc> &source, const double, const double dt) {
    const double rho = cell.density();
    for (auto d : dir::DIRECTIONS) {
        if (n(d) == 1) /* Don't integrate gravity by one-cell dimensions */
            continue;
        const double gd = g()(d);
        cell.rhoE += dt * cell.rhou(d) * gd;
        cell.rhou(d) += dt * rho * gd;
    }

    for (int i = 0; i < nc; i++)
        cell.rho[i] += dt * source.rho[i];
    cell.rhou += dt * source.rhou;
    cell.rhoE += dt * source.rhoE;
}

/* Per direction virtuals */
template<int nc>
void object<nc>::compute_inner_flux(dir::Direction dir) {
    const vec n(dir);
    int di = 0, dj = 0, dk = 0;
    dir::select(dir, di, dj, dk) = 1;
    const double hdir = h(dir);

    for (int i = di; i < nx; i++)
        for (int j = dj; j < ny; j++)
            for (int k = dk; k < nz; k++)
                flux_by(dir, i, j, k).set(
                        val(i-di, j-dj, k-dk),
                        val(i   , j   , k),
                        n, gas(), g() * hdir);
}

template<int nc>
void object<nc>::compute_outer_flux(dir::Direction dir) {
    const vec n(dir);
    const double tol = 1e-4;
    const double hdir = h(dir);
    const int ndir = dir::select(dir, nx, ny, nz);

    int ilo = 0,  jlo = 0,  klo = 0;
    int ihi = nx, jhi = ny, khi = nz;
    int di  = 0,  dj  = 0,  dk  = 0;

    dir::select(dir, ilo, jlo, klo) = 0;
    dir::select(dir, ihi, jhi, khi) = 1;

    for (int i = ilo; i < ihi; i++)
        for (int j = jlo; j < jhi; j++)
            for (int k = klo; k < khi; k++) {
                double Srefl = 1;
                flux_by(dir, i, j, k).zero();
                for (const auto &z : side(dir, dir::BEG, i, j, k)) {
                    Srefl -= z.Sfrac;
                    flux_by(dir, i, j, k).add(
                            static_cast<const object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                            val(i, j, k),
                            n, z.Sfrac, gas(), g() * hdir
                        );
                }
                if (Srefl > tol)
                    flux_by(dir, i, j, k).add_reflect(
                            nullptr,
                            val(i, j, k),
                            n, Srefl, gas(), g() * hdir);
            }

    dir::select(dir, ilo, jlo, klo) = ndir;
    dir::select(dir, ihi, jhi, khi) = ndir + 1;
    dir::select(dir, di, dj, dk) = 1;

    for (int i = ilo; i < ihi; i++)
        for (int j = jlo; j < jhi; j++)
            for (int k = klo; k < khi; k++) {
                double Srefl = 1;
                flux_by(dir, i, j, k).zero();
                for (const auto &z : side(dir, dir::END, i, j, k)) {
                    Srefl -= z.Sfrac;
                    flux_by(dir, i, j, k).add(
                            val(i - di, j - dj, k - dk),
                            static_cast<const object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                            n, z.Sfrac, gas(), g() * hdir
                        );
                }
                if (Srefl > tol)
                    flux_by(dir, i, j, k).add_reflect(
                            val(i - di, j - dj, k - dk),
                            nullptr,
                            n, Srefl, gas(), g() * hdir);
            }
}

template<int nc>
void object<nc>::integrate_by(dir::Direction dir, const double t, const double dt) {
    int di = 0, dj = 0, dk = 0;
    const double hdir = h(dir);
    dir::select(dir, di, dj, dk) = 1;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                integrate(ref(i, j, k),
                        flux_by(dir, i, j, k),
                        flux_by(dir, i+di, j+dj, k+dk),
                        dir, hdir, t, dt);
            }
}

/* Per object virtuals */
template<int nc>
double object<nc>::get_max_dt() const {
    double cmax = 0;

    double vxmax = 0;
    double vymax = 0;
    double vzmax = 0;

    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                double c = gas().sound_speed(val(i, j, k));
                vec v = val(i, j, k).velocity();
                v.x = fabs(v.x);
                v.y = fabs(v.y);
                v.z = fabs(v.z);
                if (c > cmax)
                    cmax = c;
                if (v.x > vxmax)
                    vxmax = v.x;
                if (v.y > vymax)
                    vymax = v.y;
                if (v.z > vzmax)
                    vzmax = v.z;
            }

    double dtx = h.x / (cmax + vxmax);
    double dty = h.y / (cmax + vymax);
    double dtz = h.z / (cmax + vzmax);

    return std::min(dtx, std::min(dty, dtz));
}

/* Regular methods */
template<int nc>
void object<nc>::integrate_rhs(const double t, const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                integrate_rhs(ref(i, j, k), this->source(i, j, k), t, dt);
}

template<int nc>
void object<nc>::compute_inner_fluxes() {
    compute_inner_flux(dir::X);
    compute_inner_flux(dir::Y);
    compute_inner_flux(dir::Z);
}

template<int nc>
void object<nc>::compute_outer_fluxes() {
    compute_outer_flux(dir::X);
    compute_outer_flux(dir::Y);
    compute_outer_flux(dir::Z);
}

template<int nc>
void object<nc>::integrate(const double t, const double dt) {
    integrate_by(dir::X, t, dt);
    integrate_by(dir::Y, t, dt);
    integrate_by(dir::Z, t, dt);
}

}
