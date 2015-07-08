import sys
sys.path.append('..')

import time
import math
import shutil
import os
from pyvent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40, 1.78e-5)
gas.set_component(1, 16, 1.33, 1e-5)

rho0 = 1
rho1 = .125
p0 = 100000
p1 = 10000
u0 = .75
u1 = 0
R = 8.314e3
M = 29
g = 0

# p = rho R/M T

T0 = p0 / rho0 * M / R
T1 = p1 / rho1 * M / R

class Values(Functor):
    def __call__(self, p, state):
        if p.x < -8:
            state.from_rup([30./19, 0], vec(11./30, 0, 0), 23./30, gas)
        else:
            if p.x > 10:
                state.from_rup([0.45625, 0.3], vec(0, 0, 0), 0.4, gas)
            else:
                state.from_rup([1, 0], vec(0, 0, 0), 0.4, gas)

inlet = Atm(vec(-12, -1, -1), vec(-10, 1, 1), "inlet")
pipe  = Pipe(20, 'X', vec(-10, -.5, -.5), vec(10, .5, .5), "pipe");
exit  = Room(72, 72, 72, vec(10, -4.5, -4.5), vec(19, 4.5, 4.5), "exit")

inlet.fill(Values())
pipe .fill(Values())
exit .fill(Values())

scene = [inlet, pipe, exit]

connect(inlet, pipe)
connect(pipe, exit)

solver = Solver(scene, [], .85);
solver.set_gas(gas)
solver.set_gravity(vec(0, 0, -g));

print '\n%s\n' % solver.version();

timemarks = [2 * x for x in range(100)]

prevtime = time.clock()
prevstep = 0

os.system("mkdir -p vtks/")
os.system("rm -f vtks/*")

while timemarks:
    solver.integrate()
    while timemarks and (solver.time() > timemarks[0]):
        nowtime = time.clock()
        step = solver.step();
        speed = (step - prevstep) / (nowtime - prevtime)
        lefttimesim = timemarks[-1] - solver.time(); # sim s
        timescale = speed * solver.timestep(); # sim s / rt s
        secsleft = int(lefttimesim / timescale);
        minsleft = secsleft / 60;
        secsleft = secsleft - 60 * minsleft;
        hoursleft = minsleft / 60;
        minsleft = minsleft - hoursleft * 60;
        print 'step = %8d, t = %8.2f, dt = %8.2e, speed = %6.2f steps/sec (%6.2f ms per step), %.2dh %.2dm %.2ds ETC' % \
                (step, solver.time(), solver.timestep(), speed, 1000. / speed, hoursleft, minsleft, secsleft)
        prevstep = step;
        prevtime = nowtime;
        solver.save('vtks/')
        timemarks = timemarks[1:]
