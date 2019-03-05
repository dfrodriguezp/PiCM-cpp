import numpy
from matplotlib import pyplot
import json
import os

seed = numpy.random.randint(1000000)
numpy.random.seed(seed)
os.system("mkdir -p particles_ring")
sample = open("particles_ring/sample.dat", mode="w")
N_in = 5
N_out = 100
r = 30
ring = list()
inside = list()

for i in numpy.linspace(0, 2*numpy.pi, N_out, endpoint=False):
    ring.append([r * numpy.cos(i) + r + r/10, r * numpy.sin(i) + r + r/10])

for i in range(N_in):
    theta = numpy.random.uniform(2*numpy.pi)
    random_r = numpy.random.uniform(r - r/10)
    inside.append([random_r*numpy.cos(theta) + r + r/10, random_r*numpy.sin(theta) + r + r/10])


ring = numpy.array(ring)
inside = numpy.array(inside)

for i in range(len(ring)):
    sample.write("{} {} 0.0 0.0 0.0 -1.0 0\n".format(ring[i,0], ring[i,1]))

for i in range(len(inside)):
    sample.write("{} {} 0.0 0.0 0.0 -1.0 1\n".format(inside[i,0], inside[i,1]))

sample.close()

json_data = {
    "N": N_in + N_out,
    "steps": 1000,
    "grid_size": [64, 64],
    "ss_frequency": 10,
    "dt": 0.05,
    "output": "ring_particles",
    "sample": "ring_particles/sample.dat",
    "results": ["space", "potential"],
    "seed": seed
}

with open("ring.json", "w") as data:
    json.dump(json_data, data)
