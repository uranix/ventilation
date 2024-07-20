import sys
import time
import math
import shutil
import os
from pyvent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40, 1.78e-5)
gas.set_component(1, 16, 1.20, 1e-5)

p0 = 101325
R = 8.314e3
T = 293
M0 = 29
M1 = 16
g = 0

rho0 = M0 * p0 / (R * T)
rho2 = 0.3
rho1 = rho0 - rho2 * M0 / M1

class InValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([rho0, 0.], vec(10, 0, 0), T, gas)

class Values(Functor):
    def __call__(self, p, state):
        if p.x <= 6 and p.x >= 1:
            state.from_ruT([rho1, rho2], vec(0, 0, 0), T, gas)
        else:
            state.from_ruT([rho0, 0.001 * math.sin(100 * p.x + 413 * p.y + 141 * p.z)], vec(0, 0, 0), T, gas)

F = 20
FY = 1
H = 1. / F

NL = int(sys.argv[1])

atm0   = Atm(vec(-1, -.5, -.5), vec(0, .5, .5), "atm0");
pipe1  = Pipe(10 * F - NL, 'X', vec(0, -.5, -.5), vec(10 - H * NL, .5, .5), "pipe1");
if NL > 0:
    lay2   = Room(NL, FY, F, vec(10 - H * NL, -.5, -.5), vec(10, .5, .5), "lay2");
corn3  = Room(F, FY, F, vec(10, -.5, -.5), vec(11, .5, .5), "corn3");
if NL > 0:
    lay4   = Room(F, FY, NL, vec(10, -.5, .5), vec(11, .5, .5 + H * NL), "lay4");
pipe5  = Pipe(10 * F - NL, 'Z', vec(10, -.5, .5 + H * NL), vec(11, .5, 10.5), "pipe5");
atm6   = Atm(vec(10, -.5, 10.5), vec(11, .5, 11.5), "atm6");

atm0  .fill(InValues())
pipe1 .fill(Values())
if NL > 0:
    lay2  .fill(Values())
corn3 .fill(Values())
if NL > 0:
    lay4  .fill(Values())
pipe5 .fill(Values())
atm6  .fill(Values())

connect(atm0, pipe1)
if NL > 0:
    connect(pipe1, lay2)
    connect(lay2, corn3)
    connect(corn3, lay4)
    connect(lay4, pipe5)
else:
    connect(pipe1, corn3)
    connect(corn3, pipe5)
connect(pipe5, atm6)

if NL > 0:
    scene = [atm0, pipe1, lay2, corn3, lay4, pipe5, atm6]
else:
    scene = [atm0, pipe1, corn3, pipe5, atm6]

path = [vec(0, 0, 0), vec(10.5, 0, 0), vec(10.5, 0, 10.5)]
tracer = Tracer(path, "../csv.%d/path" % NL, .05)

solver = Solver(scene, [tracer], .7);
solver.set_gas(gas)
solver.set_gravity(vec(0, 0, -g));

print('\n%s\n' % solver.version());

timemarks = [.02 * x for x in range(100)]

prevtime = time.time()
prevstep = 0

try:
    shutil.rmtree("csv.%d/" % NL)
except:
    pass

try:
    shutil.rmtree("vtks.%d/" % NL)
except:
    pass

os.mkdir("vtks.%d/" % NL)
os.mkdir("csv.%d/" % NL)

tend = timemarks[-1]
print('Simulating up to %f s' % tend)

solver.save('vtks.%d/' % NL)
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
        solver.save('vtks.%d/' % NL)
        timemarks = timemarks[1:]
