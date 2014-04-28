template<int nc>
void pipe<nc>::compute_outer_fluxes() {
	const double tol = 1e-4;
	
	if (dir == 0) {
		int j = 0;
		int k = 0;

		vec n(1, 0, 0);
		double Sfrac = 1;
		int i = 0;
		
		this->x_flux(i, j, k).zero();
		for (auto &z : this->side(0, 0, i, j, k)) {
			Sfrac -= z.Sfrac;
			this->x_flux(i, j, k).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, this->gas);
		}
		if (Sfrac > tol)
			this->x_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, this->gas);

		i = this->nx;
		Sfrac = 1;
		this->x_flux(i, j, k).zero();
		for (auto &z : this->side(0, 1, i, j, k)) {
			Sfrac -= z.Sfrac;
			this->x_flux(i, j, k).add((*this)(i-1, j, k), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, this->gas);
		}
		if (Sfrac > tol)
			this->x_flux(i, j, k).add_reflect((*this)(i-1, j, k), true, n, Sfrac, this->gas);
	}
	
	if (dir == 1) {
		int i = 0;
		int k = 0;
		vec n(0, 1, 0);
		double Sfrac = 1;
		int j = 0;
		
		this->y_flux(i, j, k).zero();
		for (auto &z : this->side(1, 0, i, j, k)) {
			Sfrac -= z.Sfrac;
			this->y_flux(i, j, k).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, this->gas);
		}
		if (Sfrac > tol)
			this->y_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, this->gas);

		j = this->ny;
		Sfrac = 1;
		this->y_flux(i, j, k).zero();
		for (auto &z : this->side(1, 1, i, j, k)) {
			Sfrac -= z.Sfrac;
			this->y_flux(i, j, k).add((*this)(i, j-1, k), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, this->gas);
		}
		if (Sfrac > tol)
			this->y_flux(i, j, k).add_reflect((*this)(i, j-1, k), true, n, Sfrac, this->gas);
	}
	if (dir == 2) {
		int i = 0;
		int j = 0;

		vec n(0, 0, 1);
		double Sfrac = 1;
		int k = 0;
		
		this->z_flux(i, j, k).zero();
		for (auto &z : this->side(2, 0, i, j, k)) {
			Sfrac -= z.Sfrac;
			this->z_flux(i, j, k).add((*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), (*this)(i, j, k), n, z.Sfrac, this->gas);
		}
		if (Sfrac > tol)
			this->z_flux(i, j, k).add_reflect((*this)(i, j, k), false, n, Sfrac, this->gas);

		k = this->nz;
		Sfrac = 1;
		this->z_flux(i, j, k).zero();
		for (auto &z : this->side(2, 1, i, j, k)) {
			Sfrac -= z.Sfrac;
			this->z_flux(i, j, k).add((*this)(i, j, k-1), (*static_cast<scene_object<nc> *>(z.other))(z.ri, z.rj, z.rk), n, z.Sfrac, this->gas);
		}
		if (Sfrac > tol)
			this->z_flux(i, j, k).add_reflect((*this)(i, j, k-1), true, n, Sfrac, this->gas);
	}
}

template<int nc>
double pipe<nc>::get_max_dt() const {
	double maxv = 0;
	double hdir = 0;
	if (dir == 0) {
		for (int i = 0; i <= this->nx; i++) {
			double v = this->x_flux(i, 0, 0).vmax;
			if (v > maxv)
				maxv = v;
		}
		hdir = this->h.x;
	}
	if (dir == 1) {
		for (int j = 0; j <= this->ny; j++) {
			double v = this->y_flux(0, j, 0).vmax;
			if (v > maxv)
				maxv = v;
		}
		hdir = this->h.y;
	}
	if (dir == 2) { 
		for (int k = 0; k <= this->nz; k++) {
			double v = this->z_flux(0, 0, k).vmax;
			if (v > maxv)
				maxv = v;
		}
		hdir = this->h.z;
	}
	return hdir / maxv;
}

template<int nc>
void pipe<nc>::integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, double dt) {
	for (int i = 0; i < nc; i++)
		cell.rho[i] -= dt * (right.fdens[i] - left.fdens[i]) / h;
	if (dir == 0) 
		cell.rhou.x -= dt * (right.fmom.x - left.fmom.x) / h;
	if (dir == 1) 
		cell.rhou.y -= dt * (right.fmom.y - left.fmom.y) / h;
	if (dir == 2) 
		cell.rhou.z -= dt * (right.fmom.z - left.fmom.z) / h;
	cell.rhoE -= dt * (right.fener - left.fener) / h;
}

template<int nc>
void pipe<nc>::integrate(const double dt) {
	for (int i = 0; i < this->nx; i++)
		for (int j = 0; j < this->ny; j++)
			for (int k = 0; k < this->nz; k++) {
				if (dir == 0) integrate((*this)(i, j, k), this->x_flux(i, j, k), this->x_flux(i+1, j, k), this->h.x, dt);
				if (dir == 1) integrate((*this)(i, j, k), this->y_flux(i, j, k), this->y_flux(i, j+1, k), this->h.y, dt);
				if (dir == 2) integrate((*this)(i, j, k), this->z_flux(i, j, k), this->z_flux(i, j, k+1), this->h.z, dt);
			}
}
