import sys

import time
import math
import shutil
import os
from pyvent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40, 1.78e-5)
gas.set_component(1, 16, 1.33, 1e-5)

p0 = 101325
R = 8.314e3
T = 293
M0 = 29
M1 = 16
g = 0
k = M0 * g / (R * T)
rCH40 = 0
rAir0 = M0 * p0 / (R * T)
rCH41 = 0.3
rAir1 = rAir0 - M0 / M1 * rCH41

F = 10
FY = 4
H = 20.

class AtmClean(Functor):
    def __call__(self, p, state):
        state.from_ruT([rAir0, rCH40], vec(1, 0, 0), T, gas)

class AtmOut(Functor):
    def __call__(self, p, state):
        state.from_ruT([rAir0, rCH40], vec(0, 0, 0), T, gas)

class Pipe(Functor):
    def __call__(self, p, state):
        if p.x > -5 or p.x < -H/2:
            state.from_ruT([rAir0, rCH40], vec(0, 0, 0), T, gas)
        else:
            state.from_ruT([rAir1, rCH41], vec(0, 0, 0), T, gas)



atm2 = Atm(vec(-1, -1, H), vec(1, 1, H+1), "out")
atm2.fill(AtmOut())
col = Room(F, FY, int(H/2*F), vec(-1, -1, 0), vec(1, 1, H), "col")
col.fill(Pipe())
turn = Room(F, FY, F, vec(-1, -1, -2), vec(1, 1, 0), "turn");
turn.fill(Pipe());
pipe = Room(int(H/4 * F), FY, F, vec(-1 - H/2, -1, -2), vec(-1, 1, 0), "pipe");
pipe.fill(Pipe());
fan = Fan('X', vec(-2 - H/2, -1, -2), vec(-1 - H/2, 1, 0), "gv", 20);
fan.fill(Pipe());
atm0 = Atm(vec(-3 - H/2, -1, -2), vec(-2 - H/2, 1, 0), "in");
atm0.fill(AtmOut());

connect(atm2, col)
connect(col, turn)
connect(turn, pipe)
connect(pipe, fan)
connect(fan, atm0)

scene = [atm2, atm0, pipe, turn, col, fan]

solver = Solver(scene, [], .25);
solver.set_gas(gas)
solver.set_gravity(vec(0, 0, -g));

timemarks = [.03 * x for x in range(1, 200)]

prevtime = time.time()
prevstep = 0

try:
    shutil.rmtree("vtks/")
except:
    pass

try:
    shutil.rmtree("csv/")
except:
    pass

os.mkdir("vtks/")
os.mkdir("csv/")

tend = timemarks[-1]
print(solver.version())
print('Simulating up to %f s' % tend)

solver.save('vtks/')
while timemarks:
    solver.integrate()
    while timemarks and (solver.time() >= timemarks[0]):
        nowtime = time.time()
        step = solver.step();
        speed = (step - prevstep) / (nowtime - prevtime)
        lefttimesim = tend - solver.time(); # sim s
        timescale = speed * solver.timestep(); # sim s / rt s
        secsleft = int(lefttimesim / timescale);
        minsleft = secsleft // 60;
        secsleft = secsleft - 60 * minsleft;
        hoursleft = minsleft // 60;
        minsleft = minsleft - hoursleft * 60;
        print('step = %8d, t = %8.2e, dt = %8.2e, speed = %6.2f steps/sec (%6.2f ms per step), %.2dh %.2dm %.2ds ETC' % \
                (step, solver.time(), solver.timestep(), speed, 1000. / speed, hoursleft, minsleft, secsleft))
        prevstep = step;
        prevtime = nowtime;
        solver.save('vtks/')
        timemarks = timemarks[1:]
