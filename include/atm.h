#ifndef __ATM_H__
#define __ATM_H__

#include "object.h"

namespace objects {

template<int nc>
struct atm : public object<nc> {
    atm(const vec &ll, const vec &ur, const std::string &id)
        : object<nc>(1, 1, 1, ll, ur, id)
    { }

    virtual void integrate_rhs(state<nc> &, const state<nc> &,
            const double, const double) override { }

    virtual void compute_outer_flux(dir::Direction) override { }
    virtual double get_max_dt() const override { return object<nc>::timestep_unconstrained; }
    virtual void integrate_by(dir::Direction, const double, const double) override { }
};

}

#endif
