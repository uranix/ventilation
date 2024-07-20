/* vim: set ft=cpp: */
%module(directors="1") pyvent
%{
#include "vec.h"
#include "state.h"
#include "object.h"
#include "room.h"
#include "pipe.h"
#include "atm.h"
#include "fan.h"
#include "pulse.h"
#include "solver.h"
%}

%include "std_string.i"
%include "std_vector.i"

%ignore operator <<(std::ostream &, const vec &);
%ignore operator *(double, const vec &);
%ignore dir::DirectionIterator;
%ignore dir::DirectionIterator::operator++();
%ignore vec::operator()(dir::Direction) const;
%ignore vec::operator=(const vec&);

%include "vec.h"

%rename(GasInfo) gasinfo;
struct gasinfo {
    gasinfo();
    %rename(set_component) set(int, double, double, double);
    void set(int nc, double molar_mass, double gamma_ratio, double viscosity);
};

%nodefaultctor box;
struct box {
};

%template(DoubleVector) std::vector<double>;

%rename(State) state;
struct state {
    void from_rue(const std::vector<double> &r, const vec &u, double eps);
    void from_rup(const std::vector<double> &r, const vec &u, double p, const gasinfo &gas);
    void from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo &gas);
};

%feature("director") functor;
%rename(Functor) functor;
struct functor {
    virtual void operator()(const vec &, state &) const = 0;
    virtual ~functor();
};

%nodefaultctor  objects::object;
%rename(Object) objects::object;
%rename(Room)   objects::room;
%rename(Pipe)   objects::pipe;
%rename(Fan)    objects::fan;
%rename(Atm)    objects::atm;
%rename(Pulse)  objects::pulse;

namespace objects {
    struct object : public box {
        void fill(const functor &f);
        void fill_sources(const functor &f);
    };
    struct room : public object {
        room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id);
    };
    struct pipe : public object {
        pipe(int n, char dir, const vec &ll, const vec &ur, const std::string &id, double friction_coeff = 0);
    };
    struct fan : public pipe {
        fan(char dir, const vec &ll, const vec &ur, const std::string &id, double Q, double friction_coeff = 0);
    };
    struct atm : public object {
        atm(const vec &ll, const vec &ur, const std::string &id);
    };
    struct pulse : public object {
        pulse(const vec &ll, const vec &ur, const std::string &id, const double drE, const double freq);
    };
}

%template(Chain)  std::vector<vec>;
%rename(Tracer) tracer;
class tracer {
public:
    tracer(const std::vector<vec> &chain, const std::string &name, double step);
};

%template(Scene) std::vector<objects::object *>;
%template(Tracers) std::vector<tracer *>;

%rename(Solver) solver;
struct solver {
    solver(
        const std::vector<objects::object *> &scene,
        const std::vector<tracer *> &tracers,
        const double cou);
    void set_gas(const gasinfo &gas);
    void set_gravity(const vec &g);
    void integrate();
    void save(const std::string &);
    std::string version() const;
    double time() const;
    double timestep() const;
    int step() const;
};

void connect(box &, box &);
