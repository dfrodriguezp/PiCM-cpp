import numpy
from matplotlib import pyplot
import click
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import json
from tqdm import tqdm
import os

@click.command()
@click.argument("file")
# @click.option("--file", default="parameters/e-e-random_Test_6.json")
def plot(file):
    with open(file) as json_data:
        d = json.load(json_data)
    phi_N = []
    # to_plot = [1800]
    to_plot = [i for i in range(0, d["steps"], 5)]  

    for t in to_plot:
        data_N = numpy.loadtxt(d["directory"] + "/phi/step" + str(t) + ".dat", unpack=True)
        phi_N.append(data_N)

    pp1 = PdfPages(d["directory"] + "/Electric_Potential.pdf")
    # os.system("mkdir -p %s" %(d["directory"] + "/videoPhi"))
    phi_only = []

    for j in range(len(phi_N)):
        phi_only.append(phi_N[j][2])
    phi_only = numpy.array(phi_only)

    for j in tqdm(range(len(phi_N))):
        x = phi_N[j][0].reshape((64, 64))
        y = phi_N[j][1].reshape((64, 64))
        phi = phi_N[j][2].reshape((64, 64))
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_title("Step %s" %to_plot[j])
        guacamayo = ax.pcolormesh(x, y, phi, cmap="jet", shading="gouraud")
        # im = ax.scatter(phi_N[j][0], phi_N[j][1], s=20, c=phi_N[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        bar = fig.colorbar(guacamayo, ax=ax)
        bar.set_label(r"$\phi(x,y)$", fontsize=20)
        guacamayo.set_clim(numpy.min(phi_only), numpy.max(phi_only))
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$", fontsize=20)
        ax.set_ylabel(r"$y$", fontsize=20)
        ax.set_xlim(0, 1.4)
        ax.set_ylim(0, 1.4)
        pyplot.tight_layout()
        # pyplot.savefig(d["directory"] + "/videoPhi/%s.jpg" %j)
        pyplot.savefig(pp1, format="pdf")

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # im = ax.scatter(phi_N[j][0], phi_N[j][1], phi_N[j][2], s=20, c=phi_N[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # bar = fig.colorbar(im, ax=ax)
        # bar.set_label(r"$\phi(x,y)$", fontsize=20)
        # im.set_clim(min(phi_N[j][2]), max(phi_N[j][2]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.legend(loc=0)
        # pyplot.tight_layout()
        # pyplot.savefig(pp1, format="pdf")
        pyplot.close()
    pp1.close()
if __name__ == "__main__":
    plot()