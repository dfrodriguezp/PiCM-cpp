import numpy
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import json
from tqdm import tqdm

def main():
    with open("../parameters.json") as data:
        params = json.load(data)
    steps = params["steps"]

    to_plot = list()
    for s in range(steps):
        if s % 100 == 0:
            to_plot.append(s)

    pp = PdfPages("phase_space.pdf")
    positions = list()
    velocities = list()

    for i in to_plot:
        pos = numpy.loadtxt("../results/positions/step_{}.dat".format(i), unpack=True)
        vel = numpy.loadtxt("../results/velocities/step_{}.dat".format(i), unpack=True)
        positions.append(pos)
        velocities.append(vel)

    N = int(params["N"] / 2)
    dr = params["dr"]
    gp = params["npoints"]
    colors = []
    for i in range(N):
        if i < N/2:
            colors.append("blue")
        else:
            colors.append("red")

    positions = numpy.array(positions)
    velocities = numpy.array(velocities)

    for d in tqdm(range(len(to_plot))):
        x = positions[d][0]
        vx = velocities[d][0]
        pyplot.figure()
        pyplot.scatter(x[::25], vx[::25], color=colors[::25], lw=0, s=5)
        pyplot.ylim(-15, 15)
        pyplot.xlim(0, dr*gp)
        pyplot.savefig(pp, format="pdf")
        pyplot.close()
    pp.close()
if __name__ == '__main__':
    main()