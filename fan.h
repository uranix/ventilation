#ifndef __FAN_H__
#define __FAN_H__

#include "pipe.h"
#include <iostream>

namespace objects {

template<int nc>
struct fan : public pipe<nc> {
	double Pmax, Qmax, S;
	vec gradPmax;
	fan(int n, int dir, const vec &ll, const vec &ur, const std::string &id, double Pmax, double Qmax)
		: pipe<nc>(n, dir, ll, ur, id), Pmax(Pmax), Qmax(Qmax)
	{
		vec P;
		P(dir) = Pmax;
		gradPmax = P / (ur - ll);
		S = this->h.x * this->h.y * this->h.z / this->h(dir);
	}
	~fan() { }
	
	virtual void integrate_rhs(state<nc> &cell, const state<nc> &source, double dt) {
		pipe<nc>::integrate_rhs(cell, source, dt);

		double s = cell.velocity()(this->dir) / Qmax;
		vec gradP = (1 - s) * gradPmax;

		cell.rhou += dt * gradP;
		cell.rhoE += dt * gradP.dot(cell.velocity());
	}
};

}

#endif
