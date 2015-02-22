#include <iostream>

#include "room.h"
#include "pipe.h"
#include "fan.h"
#include "atm.h"
#include "tracer.h"

#include "solver.h"

#include <fenv.h>

typedef objects::atm Atm;
typedef objects::pipe Pipe;
typedef objects::room Room;
typedef objects::fan Fan;
typedef objects::object Object;
typedef solver Solver;
typedef tracer Tracer;
typedef functor Functor;
typedef gasinfo GasInfo;
typedef state State;

struct AtmValues : public Functor {
    const GasInfo gas;
    AtmValues(const GasInfo &gas) : gas(gas) { }
    void operator()(const vec &, State &st) const {
        std::vector<double> r = {1.202, 0};
        st.from_ruT(r, vec(0, 0, 0), 293, gas);
    }
};

struct RoomValues : public Functor {
    const GasInfo gas;
    RoomValues(const GasInfo &gas) : gas(gas) { }
    void operator()(const vec &, State &st) const {
        std::vector<double> r = {1.23, 0};
        st.from_ruT(r, vec(0, 0, 0), 293, gas);
    }
};

struct PipeValues : public Functor {
    const GasInfo gas;
    PipeValues(const GasInfo &gas) : gas(gas) { }
    void operator()(const vec &, State &st) const {
        std::vector<double> r = {0.6863, 0.3};
        st.from_ruT(r, vec(0, 0, 0), 293, gas);
    }
};

int main() {
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);

    GasInfo gas;
    gas.set(0, 29, 1.40, 1.78e-5);
    gas.set(1, 16, 1.33, 1e-5);

    Atm  gv1   (vec(-402, -2, 200), vec(-398, 2, 204), "gv1");
    Fan  pipe2 (50, 'Z', vec(-401.5, -1.5, 1.5), vec(-398.5, 1.5, 200), "pipe2", -100000, -2.07 * 9);
    Room corn3 (1, 1, 1, vec(-401.5, -1.5, -1.5), vec(-398.5, 1.5, 1.5), "corn3");
    Pipe pipe4 (50, 'X', vec(-398.5, -1.5, -1.5), vec(-11.5, 1.5, 1.5), "pipe4");
    Room bra5  (1, 1, 1, vec(-11.5, -1.5, -1.5), vec(-8.5, 1.5, 1.5), "bra5");
    Pipe pipe6 (1,  'X', vec(-8.5, -1.5, -1.5), vec(-1.5, 1.5, 1.5), "pipe6");
    Room bra7  (1, 1, 1, vec(-1.5, -1.5, -1.5), vec(1.5, 1.5, 1.5), "bra7");
    Pipe pipe8 (50, 'X', vec(1.5, -1.5, -1.5), vec(398.5, 1.5, 1.5), "pipe8");
    Room corn9 (1, 1, 1, vec(398.5, -1.5, -1.5), vec(401.5, 1.5, 1.5), "corn9");
    Pipe pipe10(50, 'Z', vec(398.5, -1.5, 1.5), vec(401.5, 1.5, 200), "pipe10");
    Atm  atm11 (vec(398, -2, 200), vec(402, 2, 204), "atm11");
    Pipe pipe12(20, 'Y', vec(-1.5, 1.5, -1.5), vec(1.5, 198.5, 1.5), "pipe12");
    Room dend13(1, 1, 1, vec(-1.5, 198.5, -1.5), vec(1.5, 201.5, 1.5), "dend13");
    Pipe pipe14(10, 'X', vec(-9.5, 199.5, -.5), vec(-1.5, 200.5, .5), "pipe14");
    Room corn15(1, 1, 1, vec(-10.5, 199.5, -.5), vec(-9.5, 200.5, .5), "corn15");
    Fan  vmp16 (20, 'Y', vec(-10.5, 1.5, -.5), vec(-9.5, 199.5, .5), "vmp16", 100000, 12.5);

    gv1   .fill(AtmValues (gas));
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
    pipe12.fill(PipeValues(gas));
    dend13.fill(PipeValues(gas));
    pipe14.fill(RoomValues(gas));
    corn15.fill(RoomValues(gas));
    vmp16 .fill(RoomValues(gas));

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
    connect(bra7, pipe12);
    connect(pipe12, dend13);
    connect(dend13, pipe14);
    connect(pipe14, corn15);
    connect(corn15, vmp16);
    connect(bra5, vmp16);

    std::vector<Object *> scene = {
        &gv1, &pipe2, &corn3, &pipe4, &bra5, &pipe6,
        &bra7, &pipe8, &corn9, &pipe10, &atm11,
        &pipe12, &dend13, &pipe14, &corn15, &vmp16};

    Tracer tr({vec(0, 200, 0), vec(0, 0, 0), vec(400, 0, 0), vec(400, 0, 200)}, "line", 5);

    std::vector<Tracer *> tracers = {&tr};

    Solver solver(scene, tracers, .5);
    solver.set_gas(gas);
    solver.set_gravity(vec(0, 0, -9.81));

    while (solver.time() < 1000) {
        solver.integrate();
        if ((solver.step() % 2000) == 0) {
            std::cout
                << "t = " << solver.time()
                << " dt = " << solver.timestep()
                << " step = " << solver.step() << std::endl;
            solver.save("vtks/");
        }
    }
    return 0;
}
