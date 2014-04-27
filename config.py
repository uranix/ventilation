import sys
sys.path.append('.')

from vent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40)
gas.set_component(1, 16, 1.33)

room = Room(30, 30, 30, vec(0, 0, 0), vec(1, 1, 1), "room", gas);
pipe = Room(50, 1, 1, vec(-1, .4, .4), vec(0, .6, .6), "pipe", gas);

class RoomValues(Functor):
    def __call__(self, p, state):
        if (p.x + p.y + p.z < .5):
            state.from_ruT([1.5, 0.5], vec(0, 0, 0), 300, gas)
        else:
            state.from_ruT([1., 0.], vec(0, 0, 0), 300, gas)

z = RoomValues()

room.fill(z);
pipe.fill(z);
            
connect(pipe, room);

solver = Solver([room, pipe], .25);

while solver.time() < 10:
    solver.compute_fluxes()
    dt = solver.estimate_timestep()
    solver.integrate(dt)
    if ((solver.step() % 20) == 0):
        print 't = ', solver.time(), ', dt = ', dt, ', step = ', solver.step()
        solver.save()
