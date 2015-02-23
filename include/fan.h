#ifndef __FAN_H__
#define __FAN_H__

#include "../include/pipe.h"

namespace objects {

struct fan : public pipe {
    double Pmax, Qmax;
    vec gradPmax;
    fan(int n, char cdir, const vec &ll, const vec &ur, const std::string &id,
        double Pmax, double Qmax, double friction_coeff = 0
    )
        : pipe(n, cdir, ll, ur, id, friction_coeff), Pmax(Pmax), Qmax(Qmax)
    {
        vec P;
        P(this->dir) = Pmax;
        gradPmax = P / (ur - ll);
    }
    ~fan() { }

    virtual void integrate_rhs(state &cell, const state &source, const double t, const double dt) override;
};

}

#endif
