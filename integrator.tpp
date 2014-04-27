template<int nc>
void room<nc>::compute_inner_fluxes() {
	for (int i = 1; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int k = 0; k < nz; k++) {
				vec n(1, 0, 0);
				x_flux(i, j, k).zero();
				x_flux(i, j, k).add((*this)(i-1, j, k), (*this)(i, j, k), n, 1, gas);
			}
	for (int i = 0; i < nx; i++)
		for (int j = 1; j < ny; j++)
			for (int k = 0; k < nz; k++) {
				vec n(0, 1, 0);
				y_flux(i, j, k).zero();
				y_flux(i, j, k).add((*this)(i, j-1, k), (*this)(i, j, k), n, 1, gas);
			}
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int k = 1; k < nz; k++) {
				vec n(0, 0, 1);
				z_flux(i, j, k).zero();
				z_flux(i, j, k).add((*this)(i, j, k-1), (*this)(i, j, k), n, 1, gas);
			}
}

template<int nc>
void room<nc>::compute_outer_fluxes() {

	const double tol = 1e-4;

	for (int j = 0; j < ny; j++)
		for (int k = 0; k < nz; k++) {
			vec n(1, 0, 0);
			double Sfrac = 1;
			int i = 0;
			
			x_flux(i, j, k).zero();
			for (auto &z : side(0, 0, i, j, k)) {
				Sfrac -= z.Sfrac;
				x_flux(i, j, k).add((*static_cast<room<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, gas);
			}
			if (Sfrac > tol)
				x_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, gas);

			i = nx;
			Sfrac = 1;
			x_flux(i, j, k).zero();
			for (auto &z : side(0, 1, i, j, k)) {
				Sfrac -= z.Sfrac;
				x_flux(i, j, k).add((*this)(i-1, j, k), (*static_cast<room<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, gas);
			}
			if (Sfrac > tol)
				x_flux(i, j, k).add_reflect((*this)(i-1, j, k), true, n, Sfrac, gas);
		}
	for (int i = 0; i < nx; i++)
		for (int k = 0; k < nz; k++) {
			vec n(0, 1, 0);
			double Sfrac = 1;
			int j = 0;
			
			y_flux(i, j, k).zero();
			for (auto &z : side(1, 0, i, j, k)) {
				Sfrac -= z.Sfrac;
				y_flux(i, j, k).add((*static_cast<room<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, gas);
			}
			if (Sfrac > tol)
				y_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, gas);

			j = ny;
			Sfrac = 1;
			y_flux(i, j, k).zero();
			for (auto &z : side(1, 1, i, j, k)) {
				Sfrac -= z.Sfrac;
				y_flux(i, j, k).add((*this)(i, j-1, k), (*static_cast<room<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, gas);
			}
			if (Sfrac > tol)
				y_flux(i, j, k).add_reflect((*this)(i, j-1, k), true, n, Sfrac, gas);
		}
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++) {
			vec n(0, 0, 1);
			double Sfrac = 1;
			int k = 0;
			
			z_flux(i, j, k).zero();
			for (auto &z : side(2, 0, i, j, k)) {
				Sfrac -= z.Sfrac;
				z_flux(i, j, k).add((*static_cast<room<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, gas);
			}
			if (Sfrac > tol)
				z_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, gas);

			k = nz;
			Sfrac = 1;
			z_flux(i, j, k).zero();
			for (auto &z : side(2, 1, i, j, k)) {
				Sfrac -= z.Sfrac;
				z_flux(i, j, k).add((*this)(i, j, k-1), (*static_cast<room<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, gas);
			}
			if (Sfrac > tol)
				z_flux(i, j, k).add_reflect((*this)(i, j, k-1), true, n, Sfrac, gas);
		}
}

template<int nc>
double room<nc>::get_max_speed() const {
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
	return maxv;
}

template<int nc>
void room<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, double dt) {
	for (int i = 0; i < nc; i++)
		cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
	cell.rhou -= dt * (right.fmom - left.fmom) / h;
	cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void room<nc>::integrate(const double dt) {
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
			for (int k = 0; k < nz; k++) {
				integrate((*this)(i, j, k), x_flux(i, j, k), x_flux(i+1, j, k), h.x, dt);
				integrate((*this)(i, j, k), y_flux(i, j, k), y_flux(i, j+1, k), h.y, dt);
				integrate((*this)(i, j, k), z_flux(i, j, k), z_flux(i, j, k+1), h.z, dt);
			}
}
