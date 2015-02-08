#ifndef __PULSE_H__
#define __PULSE_H__

#include "object.h"

namespace objects {

template<int nc>
struct pulse : public object<nc> {
    const double drE, freq;
    double rhoE0;
    const double pi;
    pulse(const vec &ll, const vec &ur, const std::string &id, const double drE, const double freq)
        : object<nc>(1, 1, 1, ll, ur, id), drE(drE), freq(freq), pi(4 * atan(1))
    {
        rhoE0 = -1;
    }

    virtual void compute_inner_fluxes() override { }
    virtual void compute_outer_fluxes() override { }
    virtual double get_max_dt() const override { return .2 / freq; }
    virtual void integrate(const double t, const double) override {
        if (rhoE0 < 0)
            rhoE0 = this->val(0, 0, 0).rhoE;

        double A = sin(2 * pi * freq * t);
        this->ref(0, 0, 0).rhoE = rhoE0 + drE * A;
    }
    virtual void integrate_rhs(const double, const double) override { }
};

}

#endif
