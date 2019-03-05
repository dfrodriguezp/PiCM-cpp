import numpy
from matplotlib import pyplot
import json

json_data = open("ring.json", "r")
data = json.load(json_data)

step = 990
x, y = numpy.loadtxt("space/step_{}_seed_{}_.dat".format(step, data["seed"]), unpack=True)

pyplot.figure()
pyplot.scatter(x, y)
pyplot.gca().set_aspect("equal")
pyplot.show()
pyplot.close()