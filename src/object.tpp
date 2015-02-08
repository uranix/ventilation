template<int nc>
void object<nc>::compute_inner_fluxes() {
    for (int i = 1; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(1, 0, 0);
                x_flux(i, j, k).zero();
                x_flux(i, j, k).add(
                        val(i-1, j, k),
                        val(i  , j, k),
                        n, 1.0, gas(), g() * h.x);
            }
    for (int i = 0; i < nx; i++)
        for (int j = 1; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                vec n(0, 1, 0);
                y_flux(i, j, k).zero();
                y_flux(i, j, k).add(
                        val(i, j-1, k),
                        val(i, j  , k),
                        n, 1.0, gas(), g() * h.y);
            }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 1; k < nz; k++) {
                vec n(0, 0, 1);
                z_flux(i, j, k).zero();
                z_flux(i, j, k).add(
                        val(i, j, k-1),
                        val(i, j, k  ),
                        n, 1.0, gas(), g() * h.z);
            }
}

template<int nc>
void object<nc>::compute_outer_fluxes() {

    const double tol = 1e-4;

    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++) {
            vec n(1, 0, 0);
            double Sfrac = 1;
            int i = 0;

            x_flux(i, j, k).zero();
            for (auto &z : side(dir::X, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                x_flux(i, j, k).add(
                        static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                        val(i, j, k), n, z.Sfrac, gas(), g() * h.x);
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect(val(i, j, k), false, n, Sfrac, gas(), g() * h.x);

            i = nx;
            Sfrac = 1;
            x_flux(i, j, k).zero();
            for (auto &z : side(dir::X, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                x_flux(i, j, k).add(
                        val(i-1, j, k),
                        static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                        n, z.Sfrac, gas(), g() * h.x);
            }
            if (Sfrac > tol)
                x_flux(i, j, k).add_reflect(val(i-1, j, k), true, n, Sfrac, gas(), g() * h.x);
        }
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++) {
            vec n(0, 1, 0);
            double Sfrac = 1;
            int j = 0;

            y_flux(i, j, k).zero();
            for (auto &z : side(dir::Y, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add(
                        static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                        val(i, j, k), n, z.Sfrac, gas(), g() * h.y);
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect(val(i, j, k), false, n, Sfrac, gas(), g() * h.y);

            j = ny;
            Sfrac = 1;
            y_flux(i, j, k).zero();
            for (auto &z : side(dir::Y, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                y_flux(i, j, k).add(
                        val(i, j-1, k),
                        static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                        n, z.Sfrac, gas(), g() * h.y);
            }
            if (Sfrac > tol)
                y_flux(i, j, k).add_reflect(val(i, j-1, k), true, n, Sfrac, gas(), g() * h.y);
        }
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) {
            vec n(0, 0, 1);
            double Sfrac = 1;
            int k = 0;

            z_flux(i, j, k).zero();
            for (auto &z : side(dir::Z, 0, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add(
                        static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                        val(i, j, k), n, z.Sfrac, gas(), g() * h.z);
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect(val(i, j, k), false, n, Sfrac, gas(), g() * h.z);

            k = nz;
            Sfrac = 1;
            z_flux(i, j, k).zero();
            for (auto &z : side(dir::Z, 1, i, j, k)) {
                Sfrac -= z.Sfrac;
                z_flux(i, j, k).add(
                        val(i, j, k-1),
                        static_cast<object<nc> *>(z.other)->val(z.ri, z.rj, z.rk),
                        n, z.Sfrac, gas(), g() * h.z);
            }
            if (Sfrac > tol)
                z_flux(i, j, k).add_reflect(val(i, j, k-1), true, n, Sfrac, gas(), g() * h.z);
        }
}

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

template<int nc>
void object<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, dir::Direction, double h, const double, const double dt) {
    for (int i = 0; i < nc; i++)
        cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
    cell.rhou -= dt * (right.fmom - left.fmom) / h;
    cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void object<nc>::integrate(const double t, const double dt) {
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++) {
                integrate(ref(i, j, k), x_flux(i, j, k), x_flux(i+1, j, k), dir::X, h.x, t, dt);
                integrate(ref(i, j, k), y_flux(i, j, k), y_flux(i, j+1, k), dir::Y, h.y, t, dt);
                integrate(ref(i, j, k), z_flux(i, j, k), z_flux(i, j, k+1), dir::Z, h.z, t, dt);
            }
}

template<int nc>
void object<nc>::integrate_rhs(state<nc> &cell, const state<nc> &source, const double, const double dt) {
    cell.rhoE += dt * g().dot(cell.rhou);
    cell.rhou += dt * cell.density() * g();

    for (int i = 0; i < nc; i++)
        cell.rho[i] += dt * source.rho[i];
    cell.rhou += dt * source.rhou;
    cell.rhoE += dt * source.rhoE;
}

template<int nc>
void object<nc>::integrate_rhs(const double t, const double dt) {
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
void object<nc>::debug_avg() const {
    double Tavg = 0;
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                Tavg += gas().temperature(val(i, j, k));
    Tavg /= nx * ny * nz;
    std::cout << "Tavg = " << Tavg << std::endl;
}
