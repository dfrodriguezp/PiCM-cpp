import numpy
from matplotlib import pyplot
import click
from tqdm import tqdm
import json

with open("ring.json", "r") as json_data:
    data = json.load(json_data)

steps = data["steps"]
seed = data["seed"]
step = 990
phi_total = list()

px, py = numpy.loadtxt("space/step_{}_seed_{}_.dat".format(step, seed), unpack=True)
# print(len(px), len(py))
for i in tqdm(range(steps)):
    phi = numpy.loadtxt("phi/step_{}_seed_{}_.dat".format(step, seed), usecols=(2,))
    phi_total.append(phi)

x, y, phi = numpy.loadtxt("phi/step_{}_seed_{}_.dat".format(step, seed), unpack=True)
phi = phi.reshape((64, 64))
x = x.reshape((64, 64))
y = y.reshape((64, 64))

pyplot.figure()
pyplot.title(r"$\rm{step}$" + r" ${}$".format(step), fontsize=25)
d = pyplot.pcolormesh(x, y, phi, shading="gouraud")
pyplot.xlabel(r"$x \ \rm{[a.u]}$", fontsize=25)
pyplot.ylabel(r"$y \ \rm{[a.u]}$", fontsize=25)
pyplot.xlim(0, data["grid_size"][0]-1)
pyplot.ylim(0, data["grid_size"][1]-1)
bar = pyplot.colorbar(ax=pyplot.gca())
d.set_clim(numpy.min(phi_total), numpy.max(phi_total))
bar.set_label(r"$\phi \ \rm{[a.u]}$", fontsize=25)
pyplot.scatter(px, py, color="k")
pyplot.tight_layout()
pyplot.gca().set_aspect("equal")
pyplot.savefig("phi_step_{}.pdf".format(step))   
# pyplot.savefig("two_stream_{}/phi_step{}.png".format(simulation, step))   

