import sys
sys.path.append('.')

from vent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40)
gas.set_component(1, 16, 1.33)

room = Room(20, 20, 20, vec(0, 0, 0), vec(1, 1, 1), "room", gas);
pipe = Pipe(40, 0, vec(-1, .45, .45), vec(0, .55, .55), "pipe", gas);
pipe2 = Pipe(40, 1, vec(.9, 1, 0), vec(1, 1.5, .1), "pipe2", gas);
atm = Atm(vec(.8, 1.5, -.1), vec(1.1, 1.8, .2), "atm", gas);
fan = Atm(vec(-1.2, .4, .4), vec(-1, .6, .6), "fan", gas);

class FanValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([1., 0.], vec(10, 0, 0), 300, gas)

class RoomValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([0.81875, 0.1], vec(0, 0, 0), 300, gas)

class PipeValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([1., 0.], vec(0, 0, 0), 300, gas)

room.fill(RoomValues());
pipe.fill(PipeValues());
pipe2.fill(PipeValues());
atm.fill(PipeValues());
fan.fill(FanValues());
           
connect(fan, pipe);
connect(pipe, room);
connect(pipe2, room);
connect(pipe2, atm);

solver = Solver([fan, room, pipe, pipe2, atm], .25);

while solver.time() < 10:
    solver.compute_fluxes()
    dt = solver.estimate_timestep()
    solver.integrate(dt)
    if ((solver.step() % 10) == 0):
        print 't = ', solver.time(), ', dt = ', dt, ', step = ', solver.step()
        solver.save('vtks/')
