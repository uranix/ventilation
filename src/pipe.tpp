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
                static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<true>(dir),
                this->val(dir, idir).template get<false>(dir), n, z.Sfrac, this->gas(), this->g() * this->h(dir));
    }
    if (Sfrac > tol)
        this->flux_dir(dir, idir).add_reflect(
                this->val(dir, idir).template get<false>(dir), false, n, Sfrac, this->gas(), this->g() * this->h(dir));

    idir = this->n(dir);
    dir::select(dir, i, j, k) = idir;
    Sfrac = 1;
    this->flux_dir(dir, idir).zero();
    for (auto &z : this->side(dir, 1, i, j, k)) {
        Sfrac -= z.Sfrac;
        this->flux_dir(dir, idir).add(
                this->val(dir, idir - 1).template get<true>(dir),
                static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<false>(dir),
                n, z.Sfrac, this->gas(), this->g() * this->h(dir));
    }
    if (Sfrac > tol)
        this->flux_dir(dir, idir).add_reflect(
                this->val(dir, idir - 1).template get<true>(dir), true, n, Sfrac, this->gas(), this->g() * this->h(dir));

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
void pipe<nc>::integrate(sloped_state<nc> cell, const flux<nc> &left, const flux<nc> &right, dir::Direction dir, double h, const double, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.avg.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;

    cell.avg.rhou(dir) -= dt * (right.fmom(dir) - left.fmom(dir)) / h;
    cell.avg.rhoE -= dt * (right.fener - left.fener) / h;
    (void)dir;

#if SECOND_ORDER
    vec v = cell.avg.velocity();
    double p = this->gas().pressure(cell.avg);
    for (int i = 0; i < nc; i++)
        cell.slope(dir).rho[i] -= 2 * dt * (right.fdens[i] + left.fdens[i] - 2 * cell.avg.rho[i] * v(dir)) / h;
    cell.slope(dir).rhou(dir) -= 2 * dt * (right.fmom(dir) + left.fmom(dir) - 2 * (cell.avg.rhou(dir) * v(dir) + p)) / h;
    cell.slope(dir).rhoE -= 2 * dt * (right.fener + left.fener - 2 * (cell.avg.rhoE + p) * v(dir)) / h;
#endif
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
void pipe<nc>::integrate_rhs(sloped_state<nc> cell, const state<nc> &source, const double, const double dt) {
    cell.avg.rhoE += dt * this->g().dot(cell.avg.rhou);
    cell.avg.rhou(dir) += dt * cell.avg.density() * this->g()(dir);

#if SECOND_ORDER
    for (auto dir : dir::DIRECTIONS) {
        cell.slope(dir).rhoE += dt * this->g().dot(cell.slope(dir).rhou);
        cell.slope(dir).rhou(dir) += dt * cell.slope(dir).density() * this->g()(dir);
    }
#endif

    for (int i = 0; i < nc; i++)
        cell.avg.rho[i] += dt * source.rho[i];
    cell.avg.rhou += dt * source.rhou;
    cell.avg.rhoE += dt * source.rhoE;

    double Dg = 4 * surface / perimeter;
    double Re = 1 + cell.avg.rhou.norm() * Dg / this->gas().viscosity(cell.avg);
    double cf = 0.0032 + 0.221 / pow(Re, 0.237);

    cf *= friction_coeff;

    cell.avg.rhou(dir) -= dt * cell.avg.density() * cell.avg.velocity().norm2() * cf / 8;
}
