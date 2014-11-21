template<int nc>
void scene_object<nc>::compute_inner_fluxes() {
    for (int i = 1; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(1, 0, 0);
                x_flux(i, j, k).zero();
                x_flux(i, j, k).add(
                        val(i-1, j, k).template get<true >(dir::X),
                        val(i  , j, k).template get<false>(dir::X), n, 1.0, gas(), g() * h.x);
            }
    for (int i = 0; i < nx; i++)
        for (int j = 1; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(0, 1, 0);
                y_flux(i, j, k).zero();
                y_flux(i, j, k).add(
                        val(i, j-1, k).template get<true >(dir::Y),
                        val(i, j  , k).template get<false>(dir::Y), n, 1.0, gas(), g() * h.y);
            }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz; k++) {
                vec n(0, 0, 1);
                z_flux(i, j, k).zero();
                z_flux(i, j, k).add(
                        val(i, j, k-1).template get<true >(dir::Z),
                        val(i, j, k  ).template get<false>(dir::Z), n, 1.0, gas(), g() * h.z);
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
                        static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<true >(dir::X),
                        val(i, j, k).template get<false>(dir::X), n, z.Sfrac, gas(), g() * h.x);
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect(val(i, j, k).template get<false>(dir::X), false, n, Sfrac, gas(), g() * h.x);

            i = nx;
            Sfrac = 1;
            x_flux(i, j, k).zero();
            for (auto &z : side(0, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                x_flux(i, j, k).add(
                        val(i-1, j, k).template get<true >(dir::X),
                        static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<false >(dir::X),
                        n, z.Sfrac, gas(), g() * h.x);
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect(val(i-1, j, k).template get<true >(dir::X), true, n, Sfrac, gas(), g() * h.x);
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
                        static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<true >(dir::Y),
                        val(i, j, k).template get<false>(dir::Y), n, z.Sfrac, gas(), g() * h.y);
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect(val(i, j, k).template get<false>(dir::Y), false, n, Sfrac, gas(), g() * h.y);

            j = ny;
            Sfrac = 1;
            y_flux(i, j, k).zero();
            for (auto &z : side(1, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add(
                        val(i, j-1, k).template get<true >(dir::Y),
                        static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<false>(dir::Y),
                        n, z.Sfrac, gas(), g() * h.y);
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect(val(i, j-1, k).template get<true >(dir::Y), true, n, Sfrac, gas(), g() * h.y);
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
                        static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<true >(dir::Z),
                        val(i, j, k).template get<false>(dir::Z), n, z.Sfrac, gas(), g() * h.z);
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect(val(i, j, k).template get<false>(dir::Z), false, n, Sfrac, gas(), g() * h.z);

            k = nz;
            Sfrac = 1;
            z_flux(i, j, k).zero();
            for (auto &z : side(2, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add(
                        val(i, j, k-1).template get<true >(dir::Z),
                        static_cast<scene_object<nc> *>(z.other)->val(z.ri, z.rj, z.rk).template get<false>(dir::Y),
                        n, z.Sfrac, gas(), g() * h.z);
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect(val(i, j, k-1).template get<true >(dir::Y), true, n, Sfrac, gas(), g() * h.z);
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
void scene_object<nc>::integrate(sloped_state<nc> cell, const flux<nc> &left, const flux<nc> &right, dir::Direction dir, double h, const double, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.avg.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
    cell.avg.rhou -= dt * (right.fmom - left.fmom) / h;
    cell.avg.rhoE -= dt * (right.fener - left.fener) / h;
    (void)dir;

#if SECOND_ORDER
    vec v = cell.avg.velocity();
    double p = gas().pressure(cell.avg);
    for (int i = 0; i < nc; i++)
        cell.slope(dir).rho[i] -= 2 * dt * (right.fdens[i] + left.fdens[i] - 2 * cell.avg.rho[i] * v(dir)) / h;
    cell.slope(dir).rhou -= 2 * dt * (right.fmom + left.fmom - 2 * (cell.avg.rhou * v(dir) + p * dir::to_vec(dir))) / h;
    cell.slope(dir).rhoE -= 2 * dt * (right.fener + left.fener - 2 * (cell.avg.rhoE + p) * v(dir)) / h;
#endif
}

template<int nc>
void scene_object<nc>::integrate(const double t, const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                integrate(ref(i, j, k), x_flux(i, j, k), x_flux(i+1, j, k), dir::X, h.x, t, dt);
                integrate(ref(i, j, k), y_flux(i, j, k), y_flux(i, j+1, k), dir::Y, h.y, t, dt);
                integrate(ref(i, j, k), z_flux(i, j, k), z_flux(i, j, k+1), dir::Z, h.z, t, dt);
            }
}

template<int nc>
void scene_object<nc>::integrate_rhs(sloped_state<nc> cell, const state<nc> &source, const double, const double dt) {
    cell.avg.rhoE += dt * g().dot(cell.avg.rhou);
    cell.avg.rhou += dt * cell.avg.density() * g();

#if SECOND_ORDER
    for (auto dir : dir::DIRECTIONS) {
        cell.slope(dir).rhoE += dt * g().dot(cell.slope(dir).rhou);
        cell.slope(dir).rhou += dt * cell.slope(dir).density() * g();
    }
#endif

    for (int i = 0; i < nc; i++)
        cell.avg.rho[i] += dt * source.rho[i];
    cell.avg.rhou += dt * source.rhou;
    cell.avg.rhoE += dt * source.rhoE;
}

template<int nc>
void scene_object<nc>::integrate_rhs(const double t, const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                integrate_rhs(ref(i, j, k), this->source(i, j, k), t, dt);
}

double minmod(double a, double b) {
    if (a * b <= 0)
        return 0;
    if (a > 0)
        return std::min(a, b);
    return std::max(a, b);
}

double minmod3(double a, double b, double c) {
    const double theta = 1;
    return minmod(a, 0.5 * theta * minmod(b, c));
}

template<int nc>
void limit(state<nc> &slope, const state<nc> &lf, const state<nc> &ce, const state<nc> &rt) {
    for (int i = 0; i < nc; i++)
        slope.rho[i] = minmod3(slope.rho[i], ce.rho[i] - lf.rho[i], rt.rho[i] - ce.rho[i]);
    for (auto d : dir::DIRECTIONS)
        slope.rhou(d) = minmod3(slope.rhou(d), ce.rhou(d) - lf.rhou(d), rt.rhou(d) - ce.rhou(d));
    slope.rhoE = minmod3(slope.rhoE, ce.rhoE - lf.rhoE, rt.rhoE - ce.rhoE);
}

template<int nc>
void scene_object<nc>::limit_slopes() {
#if SECOND_ORDER
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            ref(0, j, k).sx.zero();
            ref(nx - 1, j, k).sx.zero();
        }
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++) {
            ref(i, 0, k).sy.zero();
            ref(i, ny - 1, k).sy.zero();
        }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
            ref(i, j, 0).sz.zero();
            ref(i, j, nz - 1).sz.zero();
        }
    for (int i = 1; i < nx - 1; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                limit(ref(i, j, k).sx, val(i-1, j, k).ce(), val(i, j, k).ce(), val(i+1, j, k).ce());
    for (int i = 0; i < nx; i++)
        for (int j = 1; j < ny - 1; j++)
            for (int k = 0; k < nz; k++)
                limit(ref(i, j, k).sy, val(i, j-1, k).ce(), val(i, j, k).ce(), val(i, j+1, k).ce());
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz - 1; k++)
                limit(ref(i, j, k).sz, val(i, j, k-1).ce(), val(i, j, k).ce(), val(i, j, k+1).ce());
#endif
}

template<int nc>
void scene_object<nc>::debug_avg() const {
    double Tavg = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                Tavg += gas().temperature(val(i, j, k).ce());
    Tavg /= nx * ny * nz;
    std::cout << "Tavg = " << Tavg << std::endl;
}
