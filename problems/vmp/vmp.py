import sys

import time
import math
import shutil
import os
from pyvent import *

gas = GasInfo()
gas.set_component(0, 29, 1.40, 1.78e-5)
gas.set_component(1, 16, 1.33, 1e-5)

F = 3

scene = []

atm0  = Atm(vec(-402, -2, 100), vec(-398, 2, 104), "atm0")
gv1   = Fan('Z', vec(-401.5, -1.5, 98), vec(-398.5, 1.5, 100), "gv1", -9 * 2.09)
pipe2 = Pipe(100 * F, 'Z', vec(-401.5, -1.5, 4.5), vec(-398.5, 1.5, 98), "pipe2")
corn3 = Room(3 * F, 3 * F, 3 * F, vec(-401.5, -1.5, 1.5), vec(-398.5, 1.5, 4.5), "corn3")
corn4 = Room(6 * F, 3 * F, 3 * F, vec(-401.5, -1.5, -1.5), vec(-395.5, 1.5, 1.5), "corn4")
pipe5 = Pipe(200 * F, 'X', vec(-395.5, -1.5, -1.5), vec(-13.5, 1.5, 1.5), "pipe5")
conn6 = Room(20 * F, 3 * F, 3 * F, vec(-13.5, -1.5, -1.5), vec(6.5, 1.5, 1.5), "conn6");
pipe7 = Pipe(200 * F, 'X', vec(6.5, -1.5, -1.5), vec(395.5, 1.5, 1.5), "pipe7")
corn8 = Room(3 * F, 3 * F, 3 * F, vec(395.5, -1.5, -1.5), vec(398.5, 1.5, 1.5), "corn8")
corn9 = Room(3 * F, 3 * F, 6 * F, vec(398.5, -1.5, -1.5), vec(401.5, 1.5, 4.5), "corn9")
pipe10= Pipe(100 * F, 'Z', vec(398.5, -1.5, 4.5), vec(401.5, 1.5, 100), "pipe10")
atm11 = Atm(vec(398, -2, 100), vec(402, 2, 104), "atm11");
vmp12 = Fan('Y', vec(-10.5, 1.5, -.5), vec(-9.5, 4.5, .5), "vmp12", 12.5)
pipe13= Pipe(100 * F, 'Y', vec(-10.5, 4.5, -.5), vec(-9.5, 195.5, .5), "pipe13")
corn14= Room(F, 3 * F, F, vec(-10.5, 195.5, -.5), vec(-9.5, 198.5, .5), "corn14")
corn15= Room(9 * F, F, F, vec(-10.5, 198.5, -.5), vec(-1.5, 199.5, .5), "corn15")
corn16= Room(3 * F, 8 * F, 3 * F, vec(-1.5, 192.5, -1.5), vec(1.5, 200.5, 1.5), "corn16")
pipe17= Pipe(100 * F, 'Y', vec(-1.5, 4.5, -1.5), vec(1.5, 192.5, 1.5), "pipe17")
bra18 = Room(3 * F, 3 * F, 3 * F, vec(-1.5, 1.5, -1.5), vec(1.5, 4.5, 1.5), "bra18")

scene.append(atm0  );
scene.append(gv1   );
scene.append(pipe2 );
scene.append(corn3 );
scene.append(corn4 );
scene.append(pipe5 );
scene.append(conn6 );
scene.append(pipe7 );
scene.append(corn8 );
scene.append(corn9 );
scene.append(pipe10);
scene.append(atm11 );
scene.append(vmp12 );
scene.append(pipe13);
scene.append(corn14);
scene.append(corn15);
scene.append(corn16);
scene.append(pipe17);
scene.append(bra18 );

p0 = 101325
R = 8.314e3
T = 293
M0 = 29
M1 = 16
g = 9.81
H = 100
k = M0 * g / (R * T)
rCH4 = 0.3

class AtmValues(Functor):
    def __call__(self, p, state):
        state.from_ruT([M0 * p0 / (R * T), 0.], vec(0, 0, 0), T, gas)

class Column(Functor):
    def __call__(self, p, state):
        h = H - p.z
        state.from_ruT([M0 * p0 / (R * T) * math.exp(k * h), 0.], vec(0, 0, 0), T, gas)

class RoomValues(Functor):
    def __call__(self, p, state):
        Z = p0 / (R * T) * math.exp(k * H)
        state.from_ruT([M0 * Z, 0], vec(0, 0, 0), T, gas)

class DeadendValues(Functor):
    def __call__(self, p, state):
        Z = p0 / (R * T) * math.exp(k * H)
        state.from_ruT([M0 * (Z - rCH4 / M1), rCH4], vec(0, 0, 0), T, gas)

atm0  .fill(AtmValues())
gv1   .fill(AtmValues())
pipe2 .fill(Column())
corn3 .fill(RoomValues())
corn4 .fill(RoomValues())
pipe5 .fill(RoomValues())
conn6 .fill(RoomValues())
pipe7 .fill(RoomValues())
corn8 .fill(RoomValues())
corn9 .fill(RoomValues())
pipe10.fill(Column())
atm11 .fill(AtmValues())
vmp12 .fill(RoomValues())
pipe13.fill(RoomValues())
corn14.fill(RoomValues())
corn15.fill(RoomValues())
corn16.fill(DeadendValues())
pipe17.fill(DeadendValues())
bra18 .fill(DeadendValues())

connect(atm0, gv1)
connect(gv1, pipe2)
connect(pipe2, corn3)
connect(corn3, corn4)
connect(corn4, pipe5)
connect(pipe5, conn6)
connect(conn6, pipe7)
connect(pipe7, corn8)
connect(corn8, corn9)
connect(corn9, pipe10)
connect(pipe10, atm11)
connect(conn6, vmp12)
connect(vmp12, pipe13)
connect(pipe13, corn14)
connect(corn14, corn15)
connect(corn15, corn16)
connect(corn16, pipe17)
connect(pipe17, bra18)
connect(bra18, conn6)

path5678 = [vec(0, 200, 0), vec(0, 0, 0), vec(400, 0, 0), vec(400, 0, 100)];
trace1 = Tracer(path5678, "../csv/path5678", .5);
path1278 = [vec(-400, 0, 100), vec(-400, 0, 0), vec(400, 0, 0), vec(400, 0, 100)];
trace2 = Tracer(path1278, "../csv/path1278", .5);

tracers = [trace1, trace2]

solver = Solver(scene, tracers, .85);
solver.set_gas(gas)
solver.set_gravity(vec(0, 0, -9.81));

timemarks = [.2 * x for x in range(1, 2510)]

prevtime = time.time()
prevstep = 0

def prepare(d):
    os.system('mkdir -p %s/; rm -f %s/*' % (d, d))

prepare('vtks')
prepare('csv')

tend = timemarks[-1]
print(solver.version())
print('Simulating up to %f s' % tend)

solver.save('vtks/')
while timemarks:
    solver.integrate()
#    solver.save('vtks/')
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
