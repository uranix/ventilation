#ifndef __ATM_H__
#define __ATM_H__

#include "scene_object.h"

namespace objects {

template<int nc>
struct atm : public scene_object<nc> {
	atm(const vec &ll, const vec &ur, const std::string &id, const gasinfo<nc> &gas)
		: scene_object<nc>(1, 1, 1, ll, ur, id, gas)
	{
	}

	virtual void compute_inner_fluxes() { }
	virtual void compute_outer_fluxes() { }
	virtual double get_max_dt() const { return 1e20; }
	virtual void integrate(const double dt) { }
};

}

#endif
