#ifndef __PIPE_H__
#define __PIPE_H__

#include "scene_object.h"

namespace objects {

template<int nc>
struct pipe : public scene_object<nc> {
	int dir;

	pipe(int n, int dir, const vec &ll, const vec &ur, const std::string &id, const gasinfo<nc> &gas)
		: scene_object<nc>(dir == 0 ? n : 1, dir == 1 ? n : 1, dir == 2 ? n : 1, ll, ur, id, gas), dir(dir)
	{
		for (int d = 0; d < 3; d++)
			for (int s = 0; s < 2; s++)
				this->closed[d][s] = dir != d;
	}

	virtual ~pipe() { }

	virtual void compute_outer_fluxes();

	virtual double get_max_dt() const;

	virtual void integrate(const double dt);
	virtual void integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, double dt);
};

#include "pipe.tpp"

}

#endif
