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

F = 3

atm0  = Atm(vec(-402, -2, 100), vec(-398, 2, 104), "atm0")
gv1   = Fan(1, 'Z', vec(-401.5, -1.5, 95), vec(-398.5, 1.5, 100), "gv1", -10000, -2.09)
pipe2 = Pipe(100, 'Z', vec(-401.5, -1.5, 1.5), vec(-398.5, 1.5, 95), "pipe2")
corn3 = Room(F, F, F, vec(-401.5, -1.5, -1.5), vec(-398.5, 1.5, 1.5), "corn3")
pipe4 = Pipe(100, 'X', vec(-398.5, -1.5, -1.5), vec(-11.5, 1.5, 1.5), "pipe4")
bra5  = Room(F, F, F, vec(-11.5, -1.5, -1.5), vec(-8.5, 1.5, 1.5), "bra5")
pipe6 = Pipe(10,  'X', vec(-8.5, -1.5, -1.5), vec(-1.5, 1.5, 1.5), "pipe6")
bra7  = Room(F, F, F, vec(-1.5, -1.5, -1.5), vec(1.5, 1.5, 1.5), "bra7")
pipe8 = Pipe(100, 'X', vec(1.5, -1.5, -1.5), vec(398.5, 1.5, 1.5), "pipe8")
corn9 = Room(F, F, F, vec(398.5, -1.5, -1.5), vec(401.5, 1.5, 1.5), "corn9")
pipe10= Pipe(100, 'Z', vec(398.5, -1.5, 1.5), vec(401.5, 1.5, 100), "pipe10")
atm11 = Atm(vec(398, -2, 100), vec(402, 2, 104), "atm11");
pipe12= Pipe(50, 'Y', vec(-1.5, 1.5, -1.5), vec(1.5, 198.5, 1.5), "pipe12")
dend13= Room(F, F, F, vec(-1.5, 198.5, -1.5), vec(1.5, 201.5, 1.5), "dend13")
pipe14= Pipe(10, 'X', vec(-9.5, 199.5, -.5), vec(-1.5, 200.5, .5), "pipe14")
corn15= Room(1, 1, 1, vec(-10.5, 199.5, -.5), vec(-9.5, 200.5, .5), "corn15")
vmp16 = Fan(50, 'Y', vec(-10.5, 1.5, -.5), vec(-9.5, 199.5, .5), "vmp16", 10000, 12.5)

p0 = 101325
R = 8.314e3
T = 293
M0 = 29
M1 = 16
g = 9.81
k = M0 * g / (R * T)

class AtmValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([M0 * p0 / (R * T), 0.], vec(0, 0, 0), T, gas)

class Column(Functor):
    def __call__(self, p, state):
        h = 100 - p.z
        state.from_ruT([M0 * p0 / (R * T) * math.exp(k * h), 0.], vec(0, 0, 0), T, gas)

class RoomValues(Functor):
    def __call__(self, p, state):
        Z = p0 / (R * T) * math.exp(k * 100)
        state.from_ruT([M0 * Z, 0], vec(0, 0, 0), 293, gas)

class PipeValues(Functor):
    def __call__(self, p, state):
        Z = p0 / (R * T) * math.exp(k * 100)
        state.from_ruT([M0 * (Z - 0.3 / M1), 0.3], vec(0, 0, 0), 293, gas)

atm0  .fill(AtmValues())
gv1   .fill(AtmValues())
pipe2 .fill(Column())
corn3 .fill(RoomValues())
pipe4 .fill(RoomValues())
bra5  .fill(RoomValues())
pipe6 .fill(RoomValues())
bra7  .fill(RoomValues())
pipe8 .fill(RoomValues())
corn9 .fill(RoomValues())
pipe10.fill(Column())
atm11 .fill(AtmValues())
pipe12.fill(PipeValues())
dend13.fill(PipeValues())
pipe14.fill(RoomValues())
corn15.fill(RoomValues())
vmp16 .fill(RoomValues())

connect(atm0, gv1)
connect(gv1, pipe2)
connect(pipe2, corn3)
connect(corn3, pipe4)
connect(pipe4, bra5)
connect(bra5, pipe6)
connect(pipe6, bra7)
connect(bra7, pipe8)
connect(pipe8, corn9)
connect(corn9, pipe10)
connect(pipe10, atm11)
connect(bra7, pipe12)
connect(pipe12, dend13)
connect(dend13, pipe14)
connect(pipe14, corn15)
connect(corn15, vmp16)
connect(bra5, vmp16)

scene = [
    atm0, gv1, pipe2, corn3, pipe4, bra5, pipe6, bra7, pipe8, corn9, pipe10,
    atm11, pipe12, dend13, pipe14, corn15, vmp16]

path5678 = [vec(0, 200, 0), vec(0, 0, 0), vec(400, 0, 0), vec(400, 0, 100)];
trace1 = Tracer(path5678, "../csv/path5678", .5);
path1278 = [vec(-400, 0, 100), vec(-400, 0, 0), vec(400, 0, 0), vec(400, 0, 100)];
trace2 = Tracer(path1278, "../csv/path1278", .5);

tracers = [trace1, trace2]

solver = Solver(scene, tracers, .25);
solver.set_gas(gas)
solver.set_gravity(vec(0, 0, -9.81));

timemarks = [.2 * x for x in range(1, 2510)]

prevtime = time.clock()
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
print solver.version()
print 'Simulating up to %f s' % tend

solver.save('vtks/')
while timemarks:
    solver.integrate()
    while timemarks and (solver.time() >= timemarks[0]):
        nowtime = time.clock()
        step = solver.step();
        speed = (step - prevstep) / (nowtime - prevtime)
        lefttimesim = tend - solver.time(); # sim s
        timescale = speed * solver.timestep(); # sim s / rt s
        secsleft = int(lefttimesim / timescale);
        minsleft = secsleft / 60;
        secsleft = secsleft - 60 * minsleft;
        hoursleft = minsleft / 60;
        minsleft = minsleft - hoursleft * 60;
        print 'step = %8d, t = %8.2e, dt = %8.2e, speed = %6.2f steps/sec (%6.2f ms per step), %.2dh %.2dm %.2ds ETC' % \
                (step, solver.time(), solver.timestep(), speed, 1000. / speed, hoursleft, minsleft, secsleft)
        prevstep = step;
        prevtime = nowtime;
        solver.save('vtks/')
        timemarks = timemarks[1:]
