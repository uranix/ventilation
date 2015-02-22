#ifndef __PIPE_H__
#define __PIPE_H__

#include "object.h"

namespace objects {

struct pipe : public object {
    dir::Direction dir;
    double surface, perimeter;
    double friction_coeff;

    pipe(int n, char cdir, const vec &ll, const vec &ur, const std::string &id, const double fc = 0)
        : object(
            dir::from_char(cdir) == dir::X ? n : 1,
            dir::from_char(cdir) == dir::Y ? n : 1,
            dir::from_char(cdir) == dir::Z ? n : 1, ll, ur, id),
        dir(dir::from_char(cdir)), friction_coeff(fc)
    {
        for (auto d : dir::DIRECTIONS)
            for (auto s : dir::SIDES)
                this->closed[d][s] = dir != d;

        surface = this->h.x * this->h.y * this->h.z / this->h(dir);
        perimeter = this->h.x + this->h.y + this->h.z - this->h(dir);
        perimeter *= 2;
    }

    virtual ~pipe() { }

    virtual void integrate(state &cell, const flux &left, const flux &right, dir::Direction dir, double h, const double t, double dt) override;
    virtual void integrate_rhs(state &cell, const state &source, const double t, const double dt) override;

    virtual void compute_outer_flux(dir::Direction) override;
    virtual void integrate_by(dir::Direction, const double t, const double dt) override;

    virtual double get_max_dt() const override;
};

}

#endif
