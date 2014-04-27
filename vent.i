%module(directors="1") vent
%{
#include "vec.h"
#include "state.h"
#include "room.h"
#include "solver.h"
%}

%include "std_string.i"
%include "std_vector.i"

%ignore operator <<(std::ostream &, const vec &);
%ignore operator *(double , const vec &);

%include "vec.h"

template<int nc>
struct gasinfo {
    gasinfo();
    %rename(set_component) set(int, double, double);
    void set(int nc, double molar_mass, double gamma_factor);
};

%nodefaultctor box;
struct box {
};

%template(DoubleVector) std::vector<double>;

template<int nc>
struct state {
    void from_rue(const std::vector<double> &r, const vec &u, double eps);
    void from_ruT(const std::vector<double> &r, const vec &u, double T, const gasinfo<nc> &gas);
};

%template(State) state<2>;

%feature("director") functor;

template<int nc>
struct functor {
    virtual void operator()(const vec &, state<nc> &) = 0;
    virtual ~functor();
};
%template(Functor) functor<2>;

template<int nc>
struct room : public box {
    room(int nx, int ny, int nz, const vec &ll, const vec &ur, const std::string &id, const gasinfo<nc> &gas);
	void fill(functor<nc> &f);
};

%template(Scene) std::vector<room<2> *>;

template<int nc>
struct solver {
	solver(const std::vector<room<nc> *> &scene, const double C);
	void compute_fluxes();
	double estimate_timestep();
	void integrate(const double dt);
	void save();
    double time() const;
    int step() const;
};

void connect(box &, box &);

%template(GasInfo) gasinfo<2>;
%template(Room) room<2>;
%template(Solver) solver<2>;
