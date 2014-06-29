template<int nc>
void scene_object<nc>::compute_inner_fluxes() {
    for (int i = 1; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(1, 0, 0);
                x_flux(i, j, k).zero();
                x_flux(i, j, k).add((*this)(i-1, j, k), (*this)(i, j, k), n, 1, gas(), g());
            }
    for (int i = 0; i < nx; i++)
        for (int j = 1; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(0, 1, 0);
                y_flux(i, j, k).zero();
                y_flux(i, j, k).add((*this)(i, j-1, k), (*this)(i, j, k), n, 1, gas(), g());
            }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz; k++) {
                vec n(0, 0, 1);
                z_flux(i, j, k).zero();
                z_flux(i, j, k).add((*this)(i, j, k-1), (*this)(i, j, k), n, 1, gas(), g());
            }
}

template<int nc>
void scene_object<nc>::compute_outer_fluxes() {

    const double tol = 1e-4;

    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            vec n(1, 0, 0);
            double Sfrac = 1;
            int i = 0;
            
            x_flux(i, j, k).zero();
            for (auto &z : side(0, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                x_flux(i, j, k).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, gas(), g());
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, gas(), g());

            i = nx;
            Sfrac = 1;
            x_flux(i, j, k).zero();
            for (auto &z : side(0, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                x_flux(i, j, k).add((*this)(i-1, j, k), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, gas(), g());
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect((*this)(i-1, j, k), true, n, Sfrac, gas(), g());
        }
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++) {
            vec n(0, 1, 0);
            double Sfrac = 1;
            int j = 0;
            
            y_flux(i, j, k).zero();
            for (auto &z : side(1, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, gas(), g());
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, gas(), g());

            j = ny;
            Sfrac = 1;
            y_flux(i, j, k).zero();
            for (auto &z : side(1, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add((*this)(i, j-1, k), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, gas(), g());
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect((*this)(i, j-1, k), true, n, Sfrac, gas(), g());
        }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
            vec n(0, 0, 1);
            double Sfrac = 1;
            int k = 0;
            
            z_flux(i, j, k).zero();
            for (auto &z : side(2, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, gas(), g());
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, gas(), g());

            k = nz;
            Sfrac = 1;
            z_flux(i, j, k).zero();
            for (auto &z : side(2, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add((*this)(i, j, k-1), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, gas(), g());
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect((*this)(i, j, k-1), true, n, Sfrac, gas(), g());
        }
}

template<int nc>
double scene_object<nc>::get_max_dt() const {
    double maxv = 0;
    for (int i = 0; i <= nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                double v = x_flux(i, j, k).vmax;
                if (v > maxv)
                    maxv = v;
            }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j <= ny; j++)
            for (int k = 0; k < nz; k++) {
                double v = y_flux(i, j, k).vmax;
                if (v > maxv)
                    maxv = v;
            }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k <= nz; k++) {
                double v = z_flux(i, j, k).vmax;
                if (v > maxv)
                    maxv = v;
            }
    return get_min_h() / maxv;
}

template<int nc>
void scene_object<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
    cell.rhou -= dt * (right.fmom - left.fmom) / h;
    cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void scene_object<nc>::integrate(const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                integrate((*this)(i, j, k), x_flux(i, j, k), x_flux(i+1, j, k), h.x, dt);
                integrate((*this)(i, j, k), y_flux(i, j, k), y_flux(i, j+1, k), h.y, dt);
                integrate((*this)(i, j, k), z_flux(i, j, k), z_flux(i, j, k+1), h.z, dt);
            }
}

template<int nc>
void scene_object<nc>::integrate_rhs(state<nc> &cell, const state<nc> &source, const double dt) {
    cell.rhoE += dt * g().dot(cell.rhou);
    cell.rhou += dt * cell.density() * g();

    for (int i = 0; i < nc; i++)
        cell.rho[i] += dt * source.rho[i];
    cell.rhou += dt * source.rhou;
    cell.rhoE += dt * source.rhoE;
}

template<int nc>
void scene_object<nc>::integrate_rhs(const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                integrate_rhs((*this)(i, j, k), this->source(i, j, k), dt);
}

template<int nc>
void scene_object<nc>::debug_avg() const {
    double Tavg = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                Tavg += gas().temperature((*this)(i, j, k));
    Tavg /= nx * ny * nz;
    std::cout << "Tavg = " << Tavg << std::endl;
}
