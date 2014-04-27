#ifndef __STATE_H__
#define __STATE_H__

#include "vec.h"
#include <vector>

template<int nc>
struct gasinfo;

template<int nc>
struct state {
	double rho[nc];
	vec rhou;
	double rhoE;

	void from_rue(const std::vector<double> &r, const vec &u, double eps) {
		double rs = 0;
		for (int i = 0; i < nc; i++) {
			rho[i] = r[i];
			rs += rho[i];
		}

		rhou = rs * u;

		rhoE = rs * (eps + 0.5 * u.norm2());
	}

	void from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo<nc> &gas);

	double density() const {
		double sum = 0;
		for (int i = 0; i < nc; i++)
			sum += rho[i];
		return sum;
	}
	
	vec velocity() const {
		double r = density();
		return rhou / r;
	}

	double specific_energy() const {
		double r = density();
		return (rhoE - 0.5 * rhou.norm2() / r) / r;
	}
};

template<int nc>
struct gasinfo {
	double Rspecific[nc]; /* J / kg / K */
	double gamma[nc];

	gasinfo() {
	}

	void set(int i, double molar_mass, double gamma_factor) {
		const double Runi = 8314.4;

		Rspecific[i] = Runi / molar_mass;
		gamma[i] = gamma_factor;
	}

	double density(const state<nc> &st) const {
		return st.density();
	}

	double specific_energy(const state<nc> &st) const {
		return st.specific_energy();
	}

	double pressure(const state<nc> &st) const {
		return (gamma_factor(st) - 1) * st.density() * st.specific_energy();
	}

	double gamma_factor(const state<nc> &st) const {
		double Cv = 0, Cp = 0;

		for (int i = 0; i < nc; i++) {
			double x = st.rho[i] * Rspecific[i] / (gamma[i] - 1);
			Cv += x;
			Cp += x * gamma[i];
		}

		return Cp / Cv;
	}

	double heat_capacity_volume(const state<nc> &st) const {
		double Cv = 0;

		for (int i = 0; i < nc; i++) {
			double x = st.rho[i] * Rspecific[i] / (gamma[i] - 1);
			Cv += x;
		}

		return Cv / st.density();
	}

	double tempertature(const state<nc> &st) const {
		return st.specific_energy() / heat_capacity_volume(st);
	}
};
	
template<int nc>
void state<nc>::from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo<nc> &gas) {
	from_rue(r, u, 0);

	double eps = gas.heat_capacity_volume(*this) * T;

	from_rue(r, u, eps);
}
	
struct avg_params {
	double density;
	vec velocity;
	double specific_energy;
	double pressure;

	double solve(const avg_params &left, const avg_params &right, const vec &norm);
};

template<int nc>
struct flux {
	double fdens[nc];
	vec fmom;
	double fener;

	double vmax;

	flux() {
		zero();
	}

	void zero() {
		for (int i = 0; i < nc; i++)
			fdens[i] = 0;
		fmom = vec(0);
		fener = 0;
		vmax = 0;
	}

	void add(const state<nc> &left, const state<nc> &right, const vec &norm, const double Sfrac, const gasinfo<nc> &gas) {
		avg_params la, ra, iface;

		la.density = left.density();
		la.velocity = left.velocity();
		la.specific_energy = left.specific_energy();
		la.pressure = gas.pressure(left);

		ra.density = right.density();
		ra.velocity = right.velocity();
		ra.specific_energy = right.specific_energy();
		ra.pressure = gas.pressure(right);

		double _vmax = iface.solve(la, ra, norm);
		if (_vmax > vmax)
			vmax = _vmax;

		double vn = iface.velocity.dot(norm);

		double theta[nc];

		if (vn > 0) {
			for (int i = 0; i < nc; i++)
				theta[i] = left.rho[i] / la.density;
		} else {
			for (int i = 0; i < nc; i++)
				theta[i] = right.rho[i] / ra.density;
		}

		for (int i = 0; i < nc; i++)
			fdens[i] += Sfrac * iface.density * theta[i] * vn;

		fmom += Sfrac * (iface.density * vn * iface.velocity + iface.pressure * norm);
		fener += Sfrac * vn * (
				iface.density * (iface.specific_energy + .5 * iface.velocity.norm2()) + iface.pressure
			);
	}

	void add_reflect(const state<nc> &inner, bool inner_is_left, const vec &norm, const double Sfrac, const gasinfo<nc> &gas) {
		state<nc> outer = inner;
		outer.rhou = inner.rhou - 2 * norm * dot(norm, inner.rhou);
		if (inner_is_left)
			add(inner, outer, norm, Sfrac, gas);
		else
			add(outer, inner, norm, Sfrac, gas);
	}
};


#endif
