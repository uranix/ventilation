#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "room.h"

#include <vector>

template<int nc>
class solver {
	std::vector<room<nc> *> scene;
	const double C;
	double t;
	int _step;
public:
	solver(const std::vector<room<nc> *> &scene, const double C) : scene(scene), C(C) 
	{
		t = 0;
		_step = 0;
	}
	void compute_fluxes() {
		for (auto p : scene)
			p->compute_inner_fluxes();
		for (auto p : scene)
			p->compute_outer_fluxes();
	}
	double estimate_timestep() {
		double dtmin = 1e20;
		for (auto p : scene) {
			double dt = p->get_min_h() / p->get_max_speed();
			if (dt < dtmin)
				dtmin = dt;
		}
		return C * dtmin;
	}
	void integrate(const double dt) {
		for (auto p : scene)
			p->integrate(dt);

		t += dt;
		_step++;
	}
	void save() {
		for (auto p : scene)
			p->save(_step);
	}
	double time() const {
		return t;
	}
	int step() const {
		return _step;
	}
};

#endif
