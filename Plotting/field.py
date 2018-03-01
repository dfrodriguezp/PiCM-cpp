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
def plot(file):
    with open(file) as json_data:
        d = json.load(json_data)
    E_field = []
    # to_plot = [1800]
    to_plot = [i for i in range(0, d["steps"], 5)]

    for t in to_plot:
        data = numpy.loadtxt("{}/Efield/step{}.dat".format(d["directory"], t), unpack=True)
        E_field.append(data)
    # os.system("mkdir -p %s" %(d["directory"] + "/videoEfield"))
    # os.system("mkdir -p %s" %(d["directory"] + "/videoEfield_x"))
    # os.system("mkdir -p %s" %(d["directory"] + "/videoEfield_y"))
    pp = PdfPages("{}/Electric_Field.pdf".format(d["directory"]))
    # pp1 = PdfPages(d["directory"] + "/E_x.pdf")
    # pp2 = PdfPages(d["directory"] + "/E_y.pdf")
    Ex_only = []
    Ey_only = []
    Etotal_only = []
    for j in range(len(E_field)):
        Ex_only.append(E_field[j][2])
        Ey_only.append(E_field[j][3])
        Etotal_only.append(numpy.sqrt(E_field[j][2]**2 + E_field[j][3]**2))

    Ex_only = numpy.array(Ex_only)
    Ey_only = numpy.array(Ey_only)
    Etotal_only = numpy.array(Etotal_only)

    for j in tqdm(range(len(E_field))):
        x = E_field[j][0].reshape((64, 64))
        y = E_field[j][1].reshape((64, 64))
        Ex = E_field[j][2].reshape((64, 64))
        Ey = E_field[j][3].reshape((64, 64))
        Ex = numpy.array(Ex)
        Ey = numpy.array(Ey)
        Etotal = numpy.sqrt(Ex**2 + Ey**2)
        # print(Ex.shape)
        # cm_subsection = numpy.linspace(0.0, 1.0, len(x))
        # cmap = pyplot.get_cmap("jet")
        # colors = [cmap(x) for x in cm_subsection]
        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        ax.set_title("Step {}".format(to_plot[j]))
        colors = ax.pcolormesh(x, y, Etotal, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$E(x,y)$", fontsize=20)
        colors.set_clim(numpy.min(Etotal_only), numpy.max(Etotal_only))
        # ax.quiver(x, y, Ex, Ey, pivot="tail")
        ax.set_xlabel(r"$x$", fontsize=20)
        ax.set_ylabel(r"$y$", fontsize=20)
        ax.set_aspect("equal")
        # pyplot.savefig(d["directory"] + "/videoEfield/%s.jpg" %j)
        pyplot.savefig(pp, format="pdf")

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111)
        # ax.set_title("Step %s" %to_plot[j])
        # colors = ax.pcolormesh(x, y, Ex, cmap="jet", shading="gouraud")
        # bar = fig.colorbar(colors, ax=ax)
        # bar.set_label(r"$E_{x}(x,y)$", fontsize=20)
        # colors.set_clim(numpy.min(Ex_only), numpy.max(Ex_only))
        # # im = ax.scatter(E_fieldx_N[j][0], E_fieldx_N[j][1], s=20, c=E_fieldx_N[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # # bar = fig.colorbar(im, ax=ax)
        # # bar.set_label(r"$E_{x}$", fontsize=20)
        # # im.set_clim(min(E_fieldx_N[j][2]), max(E_fieldx_N[j][2]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.tight_layout()
        # pyplot.savefig(d["directory"] + "/videoEfield_x/%s.jpg" %j)
        # pyplot.savefig(pp1, format="pdf")

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # im = ax.scatter(E_field[j][0], E_field[j][1], E_field[j][2], s=20, c=E_field[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # bar = fig.colorbar(im, ax=ax)
        # bar.set_label(r"$E_{x}$", fontsize=20)
        # im.set_clim(min(E_field[j][2]), max(E_field[j][2]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.legend(loc=0)
        # pyplot.tight_layout()
        # pyplot.savefig(pp1, format="pdf")

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111)
        # ax.set_title("Step %s" %to_plot[j])
        # colors = ax.pcolormesh(x, y, Ey, cmap="jet", shading="gouraud")
        # bar = fig.colorbar(colors, ax=ax)
        # bar.set_label(r"$E_{y}(x,y)$", fontsize=20)
        # colors.set_clim(numpy.min(Ey_only), numpy.max(Ey_only))
        # # im = ax.scatter(E_fieldy_N[j][0], E_fieldy_N[j][1], s=20, c=E_fieldy_N[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # # bar = fig.colorbar(im, ax=ax)
        # # bar.set_label(r"$E_{y}$", fontsize=20)
        # # im.set_clim(min(E_fieldy_N[j][2]), max(E_fieldy_N[j][2]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.tight_layout()
        # pyplot.savefig(d["directory"] + "/videoEfield_y/%s.jpg" %j)
        # pyplot.savefig(pp2, format="pdf")

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # im = ax.scatter(E_field[j][0], E_field[j][1], E_field[j][3], s=20, c=E_field[j][3], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # bar = fig.colorbar(im, ax=ax)
        # bar.set_label(r"$E_{y}$", fontsize=20)
        # im.set_clim(min(E_field[j][3]), max(E_field[j][3]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.legend(loc=0)
        # pyplot.tight_layout()
        # pyplot.savefig(pp2, format="pdf")

        pyplot.close()
    # pp1.close()
    # pp2.close()
    pp.close()
if __name__ == "__main__":
    plot()
