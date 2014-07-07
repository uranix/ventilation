template<int nc>
void scene_object<nc>::compute_inner_fluxes() {
    for (int i = 1; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(1, 0, 0);
                x_flux(i, j, k).zero();
                x_flux(i, j, k).add(
                        (*this)(i-1, j, k).template get<true >(dir::X),
                        (*this)(i  , j, k).template get<false>(dir::X), n, 1.0, gas(), g() * h.x);
            }
    for (int i = 0; i < nx; i++)
        for (int j = 1; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(0, 1, 0);
                y_flux(i, j, k).zero();
                y_flux(i, j, k).add(
                        (*this)(i, j-1, k).template get<true >(dir::Y),
                        (*this)(i, j  , k).template get<false>(dir::Y), n, 1.0, gas(), g() * h.y);
            }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz; k++) {
                vec n(0, 0, 1);
                z_flux(i, j, k).zero();
                z_flux(i, j, k).add(
                        (*this)(i, j, k-1).template get<true >(dir::Z),
                        (*this)(i, j, k  ).template get<false>(dir::Z), n, 1.0, gas(), g() * h.z);
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
                x_flux(i, j, k).add(
                        (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk).template get<true >(dir::X),
                        (*this)(i, j, k).template get<false>(dir::X), n, z.Sfrac, gas(), g() * h.x);
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect((*this)(i, j, k).template get<false>(dir::X), false, n, Sfrac, gas(), g() * h.x);

            i = nx;
            Sfrac = 1;
            x_flux(i, j, k).zero();
            for (auto &z : side(0, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                x_flux(i, j, k).add(
                        (*this)(i-1, j, k).template get<true >(dir::X),
                        (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk).template get<false >(dir::X),
                        n, z.Sfrac, gas(), g() * h.x);
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect((*this)(i-1, j, k).template get<true >(dir::X), true, n, Sfrac, gas(), g() * h.x);
        }
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++) {
            vec n(0, 1, 0);
            double Sfrac = 1;
            int j = 0;

            y_flux(i, j, k).zero();
            for (auto &z : side(1, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add(
                        (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk).template get<true >(dir::Y),
                        (*this)(i, j, k).template get<false>(dir::Y), n, z.Sfrac, gas(), g() * h.y);
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect((*this)(i, j, k).template get<false>(dir::Y), false, n, Sfrac, gas(), g() * h.y);

            j = ny;
            Sfrac = 1;
            y_flux(i, j, k).zero();
            for (auto &z : side(1, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add(
                        (*this)(i, j-1, k).template get<true >(dir::Y),
                        (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk).template get<false>(dir::Y),
                        n, z.Sfrac, gas(), g() * h.y);
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect((*this)(i, j-1, k).template get<true >(dir::Y), true, n, Sfrac, gas(), g() * h.y);
        }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
            vec n(0, 0, 1);
            double Sfrac = 1;
            int k = 0;

            z_flux(i, j, k).zero();
            for (auto &z : side(2, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add(
                        (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk).template get<true >(dir::Z),
                        (*this)(i, j, k).template get<false>(dir::Z), n, z.Sfrac, gas(), g() * h.z);
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect((*this)(i, j, k).template get<false>(dir::Z), false, n, Sfrac, gas(), g() * h.z);

            k = nz;
            Sfrac = 1;
            z_flux(i, j, k).zero();
            for (auto &z : side(2, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add(
                        (*this)(i, j, k-1).template get<true >(dir::Z),
                        (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk).template get<false>(dir::Y),
                        n, z.Sfrac, gas(), g() * h.z);
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect((*this)(i, j, k-1).template get<true >(dir::Y), true, n, Sfrac, gas(), g() * h.z);
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
void scene_object<nc>::integrate(sloped_state<nc> cell, const flux<nc> &left, const flux<nc> &right, double h, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.avg.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
    cell.avg.rhou -= dt * (right.fmom - left.fmom) / h;
    cell.avg.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void scene_object<nc>::integrate(const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                integrate(ref(i, j, k), x_flux(i, j, k), x_flux(i+1, j, k), h.x, dt);
                integrate(ref(i, j, k), y_flux(i, j, k), y_flux(i, j+1, k), h.y, dt);
                integrate(ref(i, j, k), z_flux(i, j, k), z_flux(i, j, k+1), h.z, dt);
            }
}

template<int nc>
void scene_object<nc>::integrate_rhs(sloped_state<nc> cell, const state<nc> &source, const double dt) {
    cell.avg.rhoE += dt * g().dot(cell.avg.rhou);
    cell.avg.rhou += dt * cell.avg.density() * g();

    for (int i = 0; i < nc; i++)
        cell.avg.rho[i] += dt * source.rho[i];
    cell.avg.rhou += dt * source.rhou;
    cell.avg.rhoE += dt * source.rhoE;
}

template<int nc>
void scene_object<nc>::integrate_rhs(const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                integrate_rhs(ref(i, j, k), this->source(i, j, k), dt);
}

template<int nc>
void scene_object<nc>::debug_avg() const {
    double Tavg = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                Tavg += gas().temperature((*this)(i, j, k).ce());
    Tavg /= nx * ny * nz;
    std::cout << "Tavg = " << Tavg << std::endl;
}
