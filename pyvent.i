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

#define NC 2

%include "std_string.i"
%include "std_vector.i"

%ignore operator <<(std::ostream &, const vec &);
%ignore operator *(double, const vec &);
%ignore dir::DirectionIterator;
%ignore dir::DirectionIterator::operator++();
%ignore vec::operator()(dir::Direction) const;

%include "vec.h"

template<int nc>
struct gasinfo {
    gasinfo();
    %rename(set_component) set(int, double, double, double);
    void set(int nc, double molar_mass, double gamma_factor, double viscosity);
};
%template(GasInfo) gasinfo<NC>;

%nodefaultctor box;
struct box {
};

%template(DoubleVector) std::vector<double>;

template<int nc>
struct state {
    void from_rue(const std::vector<double> &r, const vec &u, double eps);
    void from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo<nc> &gas);
};
%template(State) state<NC>;

%feature("director") functor;

template<int nc>
struct functor {
    virtual void operator()(const vec &, state<nc> &) const = 0;
    virtual ~functor();
};
%template(Functor) functor<NC>;

%nodefaultctor objects::object;
namespace objects {
    template<int nc>
    struct object : public box {
        void fill(const functor<nc> &f);
        void fill_sources(const functor<nc> &f);
        void debug_avg();
    };
    template<int nc>
    struct room : public object<nc> {
        room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id);
    };
    template<int nc>
    struct pipe : public object<nc> {
        pipe(int n, char dir, const vec &ll, const vec &ur, const std::string &id, double friction_coeff = 0);
    };
    template<int nc>
    struct fan : public pipe<nc> {
        fan(int n, char dir, const vec &ll, const vec &ur, const std::string &id, double Pmax, double Qmax, double friction_coeff = 0);
    };
    template<int nc>
    struct atm : public object<nc> {
        atm(const vec &ll, const vec &ur, const std::string &id);
    };
    template<int nc>
    struct pulse : public object<nc> {
        pulse(const vec &ll, const vec &ur, const std::string &id, const double drE, const double freq);
    };
}
%template(Object) objects::object<NC>;
%template(Room)   objects::room<NC>;
%template(Pipe)   objects::pipe<NC>;
%template(Fan)    objects::fan<NC>;
%template(Atm)    objects::atm<NC>;
%template(Pulse)  objects::pulse<NC>;
%template(Chain)  std::vector<vec>;

%rename(Tracer) tracer;
class tracer {
public:
    tracer(const std::vector<vec> &chain, const std::string &name, double step);
};

%template(Scene) std::vector<objects::object<NC> *>;
%template(Tracers) std::vector<tracer *>;

template<int nc>
struct solver {
    solver(
        const std::vector<objects::object<nc> *> &scene,
        const std::vector<tracer *> &tracers,
        const double cou);
    void set_gas(const gasinfo<nc> &gas);
    void set_gravity(const vec &g);
    void integrate();
    void save(const std::string &);
    std::string version() const;
    double time() const;
    double timestep() const;
    int step() const;
};
%template(Solver) solver<NC>;

void connect(box &, box &);