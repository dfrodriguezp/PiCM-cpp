import numpy
from matplotlib import pyplot
import click
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import json
from tqdm import tqdm
import os

def plot_Efield(file, n):
    with open(file) as json_data:
        d = json.load(json_data)
   
    E_field = []
    to_plot = [step for step in numpy.arange(0, d["Steps"], d["Steps"] / n, dtype=int)]

    for t in to_plot:
        data = numpy.loadtxt("{}/Efield/step{}.dat".format(d["Directory"], t), unpack=True)
        E_field.append(data)

    pp = PdfPages("{}/Electric_Field.pdf".format(d["Directory"]))
    pp1 = PdfPages("{}/E_x.pdf".format(d["Directory"]))
    pp2 = PdfPages("{}/E_y.pdf".format(d["Directory"]))

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
        x = E_field[j][0].reshape((d["Grid points"], d["Grid points"]))
        y = E_field[j][1].reshape((d["Grid points"], d["Grid points"]))
        Ex = E_field[j][2].reshape((d["Grid points"], d["Grid points"]))
        Ey = E_field[j][3].reshape((d["Grid points"], d["Grid points"]))
        Ex = numpy.array(Ex)
        Ey = numpy.array(Ey)
        Etotal = numpy.sqrt(Ex**2 + Ey**2)

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$t={:.1f} \ \omega_p^{}$".format(0.1*j*(d["Steps"]/n), "{-1}"), fontsize=30, y=0.99)
        colors = ax.pcolormesh(x, y, Etotal, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$E(x,y)$", fontsize=20)
        colors.set_clim(numpy.min(Etotal_only), numpy.max(Etotal_only))
        # ax.quiver(x[::10], y[::10], Ex[::10], Ey[::10], pivot="tail")
        ax.set_xlabel("$x$ [au]", fontsize=20)
        ax.set_ylabel("$y$ [au]", fontsize=20)
        ax.set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$t={:.1f} \ \omega_p^{}$".format(0.1*j*(d["Steps"]/n), "{-1}"), fontsize=30, y=0.99)
        colors = ax.pcolormesh(x, y, Ex, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label("$E_{x}(x,y)$", fontsize=20)
        colors.set_clim(numpy.min(Ex_only), numpy.max(Ex_only))
        ax.set_xlabel("$x$ [au]", fontsize=20)
        ax.set_ylabel("$y$ [au]", fontsize=20)
        ax.set_aspect("equal")
        pyplot.savefig(pp1, format="pdf")
        pyplot.close()

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

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$t={:.1f} \ \omega_p^{}$".format(0.1*j*(d["Steps"]/n), "{-1}"), fontsize=30, y=0.99)
        colors = ax.pcolormesh(x, y, Ey, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label("$E_{y}(x,y)$", fontsize=20)
        colors.set_clim(numpy.min(Ey_only), numpy.max(Ey_only))
        ax.set_xlabel("$x$ [au]", fontsize=20)
        ax.set_ylabel("$y$ [au]", fontsize=20)
        ax.set_aspect("equal")
        pyplot.savefig(pp2, format="pdf")
        pyplot.close()

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

        # pyplot.close()
    pp1.close()
    pp2.close()
    pp.close()

def plot_phi(file, n):
    with open(file) as json_data:
        d = json.load(json_data)
    
    Potential = []
    to_plot = [step for step in numpy.arange(0, d["Steps"], d["Steps"] / n, dtype=int)]  

    for t in to_plot:
        data = numpy.loadtxt("{}/phi/step{}.dat".format(d["Directory"], t), unpack=True)
        Potential.append(data)

    pp = PdfPages("{}/Electric_Potential.pdf".format(d["Directory"]))
    phi_only = []

    for j in range(len(Potential)):
        phi_only.append(Potential[j][2])
    phi_only = numpy.array(phi_only)

    for j in tqdm(range(len(Potential))):
        x = Potential[j][0].reshape((d["Grid points"], d["Grid points"]))
        y = Potential[j][1].reshape((d["Grid points"], d["Grid points"]))
        phi = Potential[j][2].reshape((d["Grid points"], d["Grid points"]))

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$t={:.1f} \ \omega_p^{}$".format(0.1*j*(d["Steps"]/n), "{-1}"), fontsize=30, y=0.99)
        colors = ax.pcolormesh(x, y, phi, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$\phi(x,y)$", fontsize=20)
        colors.set_clim(numpy.min(phi_only), numpy.max(phi_only))
        ax.set_xlabel("$x$ [au]", fontsize=20)
        ax.set_ylabel("$y$ [au]", fontsize=20)
        ax.set_ylim(0, d["Grid points"]-1)
        ax.set_xlim(0, d["Grid points"]-1)
        ax.set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # im = ax.scatter(phi[j][0], phi[j][1], phi[j][2], s=20, c=phi[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # bar = fig.colorbar(im, ax=ax)
        # bar.set_label(r"$\phi(x,y)$", fontsize=20)
        # im.set_clim(min(phi[j][2]), max(phi[j][2]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.legend(loc=0)
        # pyplot.tight_layout()
        # pyplot.savefig(pp, format="pdf")

    pp.close()   

def plot_rho(file, n):
    with open(file) as json_data:
        d = json.load(json_data)
    density = []

    to_plot = [step for step in numpy.arange(0, d["Steps"], d["Steps"]/n, dtype=int)]  

    for t in to_plot:
        data = numpy.loadtxt("{}/rho/step{}.dat".format(d["Directory"], t), unpack=True)
        density.append(data)

    pp = PdfPages("{}/chargeDensity.pdf".format(d["Directory"]))
    rho_only = []

    for j in range(len(density)):
        rho_only.append(density[j][2])
    rho_only = numpy.array(rho_only)

    for j in tqdm(range(len(density))):
        x = density[j][0].reshape((d["Grid points"], d["Grid points"]))
        y = density[j][1].reshape((d["Grid points"], d["Grid points"]))
        rho = density[j][2].reshape((d["Grid points"], d["Grid points"]))

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$t={:.1f} \ \omega_p^{}$".format(0.1*j*(d["Steps"]/n), "{-1}"), fontsize=30, y=0.99)
        colors = ax.pcolormesh(x, y, rho, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$\rho(x, y)$", fontsize=20)
        colors.set_clim(numpy.min(rho_only), numpy.max(rho_only))
        ax.set_xlabel("$x$ [au]", fontsize=20)
        ax.set_ylabel("$y$ [au]", fontsize=20)
        ax.set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

        # fig = pyplot.figure()
        # ax = fig.add_subplot(111, projection="3d")
        # im = ax.scatter(rho[j][0], rho[j][1], rho[j][2], s=20, c=rho[j][2], lw=0, marker="s", label="Step " + str(to_plot[j]))
        # bar = fig.colorbar(im, ax=ax)
        # bar.set_label(r"$\rho$", fontsize=20)
        # im.set_clim(min(rho[j][2]), max(rho[j][2]))
        # ax.set_aspect("equal")
        # ax.set_xlabel(r"$x$", fontsize=20)
        # ax.set_ylabel(r"$y$", fontsize=20)
        # ax.set_xlim(0, 1.4)
        # ax.set_ylim(0, 1.4)
        # pyplot.legend(loc=0)
        # pyplot.tight_layout()
        # pyplot.savefig(pp, format="pdf")
        
    pp.close()    

@click.command()
@click.option("--name", type=str)
@click.option("--n", type=int)
@click.argument("file")
def main(file, name, n):
    if name == "Efield":
        plot_Efield(file, n)
    elif name == "phi":
        plot_phi(file, n)
    elif name == "rho":
        plot_rho(file, n)

if __name__ == '__main__':
    main()