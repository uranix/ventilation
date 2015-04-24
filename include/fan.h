#ifndef __FAN_H__
#define __FAN_H__

#include "../include/pipe.h"

namespace objects {

struct fan : public pipe {

    double U;
    flux FL, FR;

    fan(char cdir, const vec &ll, const vec &ur, const std::string &id,
        double Q, double friction_coeff = 0
    )
        : pipe(2, cdir, ll, ur, id, friction_coeff)
    {
        U = Q / this->surface;
    }
    ~fan() { }

    virtual void compute_special_flux(dir::Direction dir, const double dt_h) override;
    virtual void integrate(state &cell, const flux &left, const flux &right, dir::Direction, double h, const double, const double dt) override;

};

}

#endif
