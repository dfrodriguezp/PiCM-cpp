import numpy
from matplotlib import pyplot
import click
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import json
from tqdm import tqdm

@click.command()
@click.argument("file")
def plot(file):
    with open(file) as json_data:
        d = json.load(json_data)
    rho = []

    to_plot = [i for i in range(0, d["steps"], 200)]  

    for t in to_plot:
        data_N = numpy.loadtxt(d["directory"] + "/rho/step" + str(t) + ".dat", unpack=True)
        rho.append(data_N)

    pp = PdfPages(d["directory"] + "/rho.pdf")

    for j in tqdm(range(len(rho))):
        x = rho[j][0].reshape((51, 51))
        y = rho[j][1].reshape((51, 51))
        data = rho[j][2].reshape((51, 51))
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        guacamayo = ax.pcolormesh(x, y, data, cmap="jet", shading="gouraud", label="Step %s" %(to_plot[j]))
        bar = fig.colorbar(guacamayo, ax=ax)
        bar.set_label(r"$\rho$", fontsize=20)
        guacamayo.set_clim(min(rho[j][2]), max(rho[j][2]))
        # im = ax.scatter(rho[j][0], rho[j][1], s=20, c=rho[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # bar = fig.colorbar(im, ax=ax)
        # bar.set_label(r"$\rho$", fontsize=20)
        # im.set_clim(min(rho[j][2]), max(rho[j][2]))
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$", fontsize=20)
        ax.set_ylabel(r"$y$", fontsize=20)
        ax.set_xlim(0, 1.4)
        ax.set_ylim(0, 1.4)
        pyplot.legend(loc=0)
        pyplot.tight_layout()
        pyplot.savefig(pp, format="pdf")

        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection="3d")
        im = ax.scatter(rho[j][0], rho[j][1], rho[j][2], s=20, c=rho[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        bar = fig.colorbar(im, ax=ax)
        bar.set_label(r"$\rho$", fontsize=20)
        im.set_clim(min(rho[j][2]), max(rho[j][2]))
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$", fontsize=20)
        ax.set_ylabel(r"$y$", fontsize=20)
        ax.set_xlim(0, 1.4)
        ax.set_ylim(0, 1.4)
        pyplot.legend(loc=0)
        pyplot.tight_layout()
        pyplot.savefig(pp, format="pdf")
        
        pyplot.close()
    pp.close()

if __name__ == "__main__":
    plot()