#ifndef __PIPE_H__
#define __PIPE_H__

#include "scene_object.h"

namespace objects {

inline int cdir_to_dir(char cdir) {
    if (cdir == 'x' || cdir == 'X')
        return 0;
    if (cdir == 'y' || cdir == 'Y')
        return 1;
    return 2;
}

template<int nc>
struct pipe : public scene_object<nc> {
    int dir;
    double surface, perimeter;
    double friction_coeff;

    pipe(int n, char cdir, const vec &ll, const vec &ur, const std::string &id, double friction_coeff = 0)
        : scene_object<nc>(
            cdir_to_dir(cdir) == 0 ? n : 1,
            cdir_to_dir(cdir) == 1 ? n : 1,
            cdir_to_dir(cdir) == 2 ? n : 1, ll, ur, id), 
        dir(cdir_to_dir(cdir)), friction_coeff(friction_coeff)
    {
        for (int d = 0; d < 3; d++)
            for (int s = 0; s < 2; s++)
                this->closed[d][s] = dir != d;

        surface = this->h.x * this->h.y * this->h.z / this->h(dir);
        perimeter = this->h.x + this->h.y + this->h.z - this->h(dir);
        perimeter *= 2;
    }

    virtual ~pipe() { }

    virtual void compute_outer_fluxes();

    virtual double get_max_dt() const;

    virtual void integrate(const double dt);
    virtual void integrate(state<nc> &cell, const flux<nc> &left, const flux<nc> &right, double h, double dt);
    virtual void integrate_rhs(state<nc> &cell, const state<nc> &source, const double dt);
};

#include "../src/pipe.tpp"

}

#endif
