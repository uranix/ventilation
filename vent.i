%module(directors="1") vent
%{
#include "vec.h"
#include "state.h"
#include "scene_object.h"
#include "room.h"
#include "pipe.h"
#include "atm.h"
#include "fan.h"
#include "solver.h"
%}

#define NC 2

%include "std_string.i"
%include "std_vector.i"

%ignore operator <<(std::ostream &, const vec &);
%ignore operator *(double, const vec &);
%ignore DirectionIterator;
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

%nodefaultctor objects::scene_object;
namespace objects {
    template<int nc>
    struct scene_object : public box {
        void fill(const functor<nc> &f);
        void fill_sources(const functor<nc> &f);
        void debug_avg();
    };
    template<int nc>
    struct room : public scene_object<nc> {
        room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id);
    };
    template<int nc>
    struct pipe : public scene_object<nc> {
        pipe(int n, char dir, const vec &ll, const vec &ur, const std::string &id, double friction_coeff = 0);
    };
    template<int nc>
    struct fan : public pipe<nc> {
        fan(int n, char dir, const vec &ll, const vec &ur, const std::string &id, double Pmax, double Qmax, double friction_coeff = 0);
    };
    template<int nc>
    struct atm : public scene_object<nc> {
        atm(const vec &ll, const vec &ur, const std::string &id);
    };
}
%template(SceneObject) objects::scene_object<NC>;
%template(Room) objects::room<NC>;
%template(Pipe) objects::pipe<NC>;
%template(Fan) objects::fan<NC>;
%template(Atm) objects::atm<NC>;
%template(Scene) std::vector<objects::scene_object<NC> *>;

template<int nc>
struct solver {
    solver(const std::vector<objects::scene_object<nc> *> &scene, const double C);
    void set_gas(const gasinfo<nc> &gas);
    void set_gravity(const vec &g);
    void integrate();
    void save(const std::string &);
    double time() const;
    double timestep() const;
    int step() const;
};
%template(Solver) solver<NC>;

void connect(box &, box &);
