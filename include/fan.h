#ifndef __FAN_H__
#define __FAN_H__

#include "pipe.h"
#include <iostream>

namespace objects {

template<int nc>
struct fan : public pipe<nc> {
    double Pmax, Qmax;
    vec gradPmax;
    fan(int n, char cdir, const vec &ll, const vec &ur, const std::string &id,
        double Pmax, double Qmax, double friction_coeff = 0
    )
        : pipe<nc>(n, cdir, ll, ur, id, friction_coeff), Pmax(Pmax), Qmax(Qmax)
    {
        vec P;
        P(this->dir) = Pmax;
        gradPmax = P / (ur - ll);
    }
    ~fan() { }

    virtual void integrate_rhs(sloped_state<nc> cell, const state<nc> &source, double dt) {
        pipe<nc>::integrate_rhs(cell, source, dt);

        double s = this->surface * cell.avg.velocity()(this->dir) / Qmax;
        vec gradP = (1 - s) * gradPmax;

        cell.avg.rhou += dt * gradP;
        cell.avg.rhoE += dt * gradP.dot(cell.avg.velocity());
    }
};

}

#endif
