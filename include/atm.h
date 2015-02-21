#ifndef __ATM_H__
#define __ATM_H__

#include "object.h"

namespace objects {

struct atm : public object {
    atm(const vec &ll, const vec &ur, const std::string &id)
        : object(1, 1, 1, ll, ur, id)
    { }

    virtual void integrate_rhs(state &, const state &,
            const double, const double) override { }

    virtual void compute_outer_flux(dir::Direction) override { }
    virtual double get_max_dt() const override { return object::timestep_unconstrained; }
    virtual void integrate_by(dir::Direction, const double, const double) override { }
};

}

#endif
