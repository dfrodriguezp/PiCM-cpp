import numpy
import random
from itertools import product

# Electron charge and mass
q = -1.0
m = 1.0
qm = 1.0

class Particle(object):
    def __init__(self, pos, vel, n, qm, move):
        self.pos = pos
        self.vel = vel
        self.qm = qm
        self.q = (1/qm) * (1/n)
        self.move = move

def two_stream(N, L, dr, vt, vd, s):
    S1 = numpy.random.seed(s)
    S2 = random.seed(s)
    parts = list()
    n = N / (L*L)
    margin = dr / 10

    P = numpy.linspace(margin, L - margin, int(numpy.sqrt(N)))
    pos = numpy.array(list(product(P, P)))
    # for i in range(N):
    #     print(pos[i][0], pos[i][1])
        
    vel0 = numpy.array([0.0, 0.0])
    indexes = list(range(N))
    # random.shuffle(indexes)
    left = [indexes.pop() for _ in range(int(N/4))]
    right = [indexes.pop() for _ in range(int(N/4))]

    for i in indexes:
        parts.append(Particle(pos[i], vel0, n, qm, False))
    for i in left:
        # vel1 = numpy.array([numpy.random.normal(vd, vt), 0.0])
        vel1 = numpy.array([0.02, 0.0])
        parts.append(Particle(pos[i], vel1, n, -qm, True))
    for i in right:
        # vel2 = numpy.array([numpy.random.normal(-vd, vt), 0.0])
        vel2 = numpy.array([-0.02, 0.0])
        parts.append(Particle(pos[i], vel2, n, -qm, True))

    return parts