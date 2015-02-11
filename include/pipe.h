#ifndef __PIPE_H__
#define __PIPE_H__

#include "object.h"

namespace objects {

template<int nc>
struct pipe : public object<nc> {
    dir::Direction dir;
    double surface, perimeter;
    double friction_coeff;

    pipe(int n, char cdir, const vec &ll, const vec &ur, const std::string &id, double friction_coeff = 0)
        : object<nc>(
            dir::from_char(cdir) == dir::X ? n : 1,
            dir::from_char(cdir) == dir::Y ? n : 1,
            dir::from_char(cdir) == dir::Z ? n : 1, ll, ur, id),
        dir(dir::from_char(cdir)), friction_coeff(friction_coeff)
    {
        for (int d = 0; d < 3; d++)
            for (int s = 0; s < 2; s++)
                this->closed[d][s] = dir != d;

        surface = this->h.x * this->h.y * this->h.z / this->h(dir);
        perimeter = this->h.x + this->h.y + this->h.z - this->h(dir);
        perimeter *= 2;
    }

    virtual ~pipe() { }

    virtual void compute_outer_fluxes() override;

    virtual double get_max_dt() const override;

    virtual void integrate(const double t, const double dt) override ;
    virtual void integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, dir::Direction dir, double h, const double t, double dt) override;
    virtual void integrate_rhs(state<nc> &cell, const state<nc> &source, const double t, const double dt) override;
};

}

#endif
