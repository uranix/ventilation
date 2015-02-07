#ifndef __ATM_H__
#define __ATM_H__

#include "scene_object.h"

namespace objects {

template<int nc>
struct atm : public scene_object<nc> {
    atm(const vec &ll, const vec &ur, const std::string &id)
        : scene_object<nc>(1, 1, 1, ll, ur, id)
    {
    }

    virtual void compute_inner_fluxes() override { }
    virtual void compute_outer_fluxes() override { }
    virtual double get_max_dt() const override { return 1e20; }
    virtual void integrate(const double, const double) override { }
    virtual void integrate_rhs(const double, const double) override { }
};

}

#endif
