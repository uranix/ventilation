#include <iostream>

#include "room.h"
#include "pipe.h"
#include "fan.h"
#include "atm.h"

#include "solver.h"

#include <fenv.h>

struct AtmValues : public functor<2> {
    const gasinfo<2> gas;
    AtmValues(const gasinfo<2> &gas) : gas(gas) { }
    void operator()(const vec &p, state<2> &st) const {
        std::vector<double> r = {1.202, 0};
        st.from_ruT(r, vec(0, 0, 0), 293, gas);
    }
};

struct GvValues : public functor<2> {
    const gasinfo<2> gas;
    GvValues(const gasinfo<2> &gas) : gas(gas) { }
    void operator()(const vec &p, state<2> &st) const {
        std::vector<double> r = {1.202, 0};
        st.from_ruT(r, vec(0, 0, -2.07), 293, gas);
    }
};

struct RoomValues : public functor<2> {
    const gasinfo<2> gas;
    RoomValues(const gasinfo<2> &gas) : gas(gas) { }
    void operator()(const vec &p, state<2> &st) const {
        std::vector<double> r = {1.274, 0};
        st.from_ruT(r, vec(0, 0, 0), 283, gas);
    }
};

struct PipeValues : public functor<2> {
    const gasinfo<2> gas;
    PipeValues(const gasinfo<2> &gas) : gas(gas) { }
    void operator()(const vec &p, state<2> &st) const {
        std::vector<double> r = {0.73, 0.3};
        st.from_ruT(r, vec(0, 0, 0), 283, gas);
    }
};

typedef objects::atm<2> Atm;
typedef objects::pipe<2> Pipe;
typedef objects::room<2> Room;
typedef objects::fan<2> Fan;
typedef objects::scene_object<2> Object;
typedef solver<2> Solver;

int main() {
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

    gasinfo<2> gas;
    gas.set(0, 29, 1.40, 1.78e-5);
    gas.set(1, 16, 1.33, 1e-5);

    Atm  gv1   (vec(-402, -2, 200), vec(-398, 2, 204), "gv1");
    Pipe pipe2 (50, 'Z', vec(-401.5, -1.5, 1.5), vec(-398.5, 1.5, 200), "pipe2");
    Room corn3 (1, 1, 1, vec(-401.5, -1.5, -1.5), vec(-398.5, 1.5, 1.5), "corn3");
    Pipe pipe4 (10, 'X', vec(-398.5, -1.5, -1.5), vec(-11.5, 1.5, 1.5), "pipe4");
    Room bra5  (1, 1, 1, vec(-11.5, -1.5, -1.5), vec(-8.5, 1.5, 1.5), "bra5");
    Pipe pipe6 (1,  'X', vec(-8.5, -1.5, -1.5), vec(-1.5, 1.5, 1.5), "pipe6");
    Room bra7  (1, 1, 1, vec(-1.5, -1.5, -1.5), vec(1.5, 1.5, 1.5), "bra7");
    Pipe pipe8 (10, 'X', vec(1.5, -1.5, -1.5), vec(398.5, 1.5, 1.5), "pipe8");
    Room corn9 (1, 1, 1, vec(398.5, -1.5, -1.5), vec(401.5, 1.5, 1.5), "corn9");
    Pipe pipe10(50, 'Z', vec(398.5, -1.5, 1.5), vec(401.5, 1.5, 200), "pipe10");
    Atm  atm11 (vec(398, -2, 200), vec(402, 2, 204), "atm11");

    gv1   .fill(GvValues  (gas));
    pipe2 .fill(RoomValues(gas));
    corn3 .fill(RoomValues(gas));
    pipe4 .fill(RoomValues(gas));
    bra5  .fill(RoomValues(gas));
    pipe6 .fill(RoomValues(gas));
    bra7  .fill(RoomValues(gas));
    pipe8 .fill(RoomValues(gas));
    corn9 .fill(RoomValues(gas));
    pipe10.fill(RoomValues(gas));
    atm11 .fill(AtmValues (gas));

    connect(gv1, pipe2);
    connect(pipe2, corn3);
    connect(corn3, pipe4);
    connect(pipe4, bra5);
    connect(bra5, pipe6);
    connect(pipe6, bra7);
    connect(bra7, pipe8);
    connect(pipe8, corn9);
    connect(corn9, pipe10);
    connect(pipe10, atm11);

    std::vector<Object *> scene = {
        &gv1, &pipe2, &corn3, &pipe4, &bra5,
        &pipe6, &bra7, &pipe8, &corn9, &pipe10,
        &atm11};

    Solver solver(scene, .25);
    solver.set_gas(gas);
    solver.set_gravity(vec(0, 0, -9.81));

    while (solver.time() < 1000) {
        solver.compute_fluxes();
        double dt = solver.estimate_timestep();
        solver.integrate(dt);
        if ((solver.step() % 500) == 0) {
            std::cout << "t = " << solver.time() << " dt = "
                << dt << " step = " << solver.step() << std::endl;
            solver.save("vtks/");
        }
    }
    return 0;
}
