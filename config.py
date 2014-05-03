import sys
sys.path.append('.')

from vent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40, 1.78e-5)
gas.set_component(1, 16, 1.33, 1e-5)

Room.set_gas(gas);

room = Room(20, 20, 20, vec(0, 0, 0), vec(1, 1, 1), "room");
fan = Fan(10, 0, vec(-1, .45, .45), vec(0, .55, .55), "fan", 50000, 10);
pipe2 = Pipe(10, 1, vec(.9, 1, 0), vec(1, 1.5, .1), "pipe2");
atm = Atm(vec(.8, 1.5, -.1), vec(1.1, 1.8, .2), "atm");
atm2 = Atm(vec(-1.2, .4, .4), vec(-1, .6, .6), "atm2");

class AtmValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([1., 0.], vec(0, 0, 0), 300, gas)

class RoomValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([0.81875, 0.1], vec(0, 0, 0), 300, gas)

class PipeValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([1., 0.], vec(0, 0, 0), 300, gas)

room.fill(RoomValues());
pipe2.fill(PipeValues());
atm.fill (AtmValues());
atm2.fill(AtmValues());
fan.fill(PipeValues());
           
connect(fan, atm2);
connect(fan, room);
connect(pipe2, room);
connect(pipe2, atm);

solver = Solver([fan, room, pipe2, atm, atm2], .25);

while solver.time() < 10:
    solver.compute_fluxes()
    dt = solver.estimate_timestep()
    solver.integrate(dt)
    if ((solver.step() % 50) == 0):
        print 't =', solver.time(), ' dt =', dt, ' step =', solver.step()
        solver.save('vtks/')
