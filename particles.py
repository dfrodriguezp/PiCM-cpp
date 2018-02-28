import numpy
import random
from itertools import product

# Electron charge and mass
q = -1.0
m = 1.0

class Particle(object):
    def __init__(self, pos, vel, m, n, move):
        self.pos = pos
        self.vel = vel
        self.q = 1/n
        self.m = m
        self.qm = n
        self.move = move

def maxwell_distribution(n, vt, vd, vx):
    return (n / (numpy.sqrt(2 * numpy.pi)) * vt) * numpy.exp(-((vx - vd)**2) / (2 * vt**2))

def two_stream(N, L, dr, vt, vd, s):
    S1 = numpy.random.seed(s)
    S2 = random.seed(s)
    parts = list()
    n = N / (L**2)
    Ne_Ni = 2
    Ne = N / Ne_Ni
    margin = dr / 10

    P = numpy.linspace(0 + margin, L - margin, int(numpy.sqrt(N)))
    pos = numpy.array(list(product(P, P)))

    vx = numpy.linspace(-3*dr, 3*dr, 1000)
    f1 = maxwell_distribution(n, vt, vd, vx)
    f2 = maxwell_distribution(n, vt, vd, -vx)
    p1 = f1 / numpy.sum(f1)
    p2 = f2 / numpy.sum(f2)
    vel0 = numpy.array([0.0, 0.0])

    indexes = list(range(N))
    random.shuffle(indexes)
    left = [indexes.pop() for _ in range(int(Ne/2))]
    right = [indexes.pop() for _ in range(int(Ne/2))]

    for i in indexes:
        parts.append(Particle(pos[i], vel0, m, n, False))
    for i in left:
        vel1 = numpy.array([numpy.random.choice(vx, p=p1), 0.0])
        parts.append(Particle(pos[i], vel1, m, -n, True))
    for i in right:
        vel2 = numpy.array([numpy.random.choice(vx, p=p2), 0.0])
        parts.append(Particle(pos[i], vel2, m, -n, True))

    return parts
