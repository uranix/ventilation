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
        this->flux_dir(dir, idir).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk),
                (*this)(dir, idir), n, z.Sfrac, this->gas(), this->g() * this->h(dir));
    }
    if (Sfrac > tol)
        this->flux_dir(dir, idir).add_reflect((*this)(dir, idir), false, n, Sfrac, this->gas(), this->g() * this->h(dir));

    idir = this->n(dir);
    dir::select(dir, i, j, k) = idir;
    Sfrac = 1;
    this->flux_dir(dir, idir).zero();
    for (auto &z : this->side(dir, 1, i, j, k)) {
        Sfrac -= z.Sfrac;
        this->flux_dir(dir, idir).add((*this)(dir, idir - 1), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk),
                n, z.Sfrac, this->gas(), this->g() * this->h(dir));
    }
    if (Sfrac > tol)
        this->flux_dir(dir, idir).add_reflect((*this)(dir, idir - 1), true, n, Sfrac, this->gas(), this->g() * this->h(dir));

}

template<int nc>
double pipe<nc>::get_max_dt() const {
    double maxv = 0;
    double hdir = this->h(dir);

    for (int i = 0; i <= this->n(dir); i++) {
        double v = this->flux_dir(dir, i).vmax;
        if (v > maxv)
            maxv = v;
    }

    return hdir / maxv;
}

template<int nc>
void pipe<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;

    cell.rhou(dir) -= dt * (right.fmom(dir) - left.fmom(dir)) / h;

    cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void pipe<nc>::integrate(const double dt) {
    for (int i = 0; i < this->nx; i++)
        for (int j = 0; j < this->ny; j++)
            for (int k = 0; k < this->nz; k++) {
                if (dir == dir::X) integrate((*this)(i, j, k), this->x_flux(i, j, k), this->x_flux(i+1, j, k), this->h.x, dt);
                if (dir == dir::Y) integrate((*this)(i, j, k), this->y_flux(i, j, k), this->y_flux(i, j+1, k), this->h.y, dt);
                if (dir == dir::Z) integrate((*this)(i, j, k), this->z_flux(i, j, k), this->z_flux(i, j, k+1), this->h.z, dt);
            }
}

template<int nc>
void pipe<nc>::integrate_rhs(state<nc> &cell, const state<nc> &source, const double dt) {
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
