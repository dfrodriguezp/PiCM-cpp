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
    to_plot = [step for step in numpy.arange(0, d["steps"], d["steps"] / n, dtype=int)]

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
        x = E_field[j][0].reshape((d["gp"], d["gp"]))
        y = E_field[j][1].reshape((d["gp"], d["gp"]))
        Ex = E_field[j][2].reshape((d["gp"], d["gp"]))
        Ey = E_field[j][3].reshape((d["gp"], d["gp"]))
        Ex = numpy.array(Ex)
        Ey = numpy.array(Ey)
        Etotal = numpy.sqrt(Ex**2 + Ey**2)

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", 0.1*j*(d["steps"]/n)), fontsize=25, y=0.99)
        colors = ax.pcolormesh(x, y, Etotal, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$E \left(x,y \right) \left(e / m_{e} \omega_p v_{th} \right)$", fontsize=20)
        colors.set_clim(numpy.min(Etotal_only), numpy.max(Etotal_only))
        # ax.quiver(x[::10], y[::10], Ex[::10], Ey[::10], pivot="tail")
        ax.set_xlabel(r"$x / \lambda_{D}$", fontsize=20)
        ax.set_ylabel(r"$y / \lambda_{D}$", fontsize=20)
        ax.set_ylim(0, d["gp"]-1)
        ax.set_xlim(0, d["gp"]-1)
        ax.set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", 0.1*j*(d["steps"]/n)), fontsize=25, y=0.99)
        colors = ax.pcolormesh(x, y, Ex, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$E_{x} \left(x,y \right) \left(e / m_{e} \omega_p v_{th} \right)$", fontsize=20)
        colors.set_clim(numpy.min(Ex_only), numpy.max(Ex_only))
        ax.set_xlabel(r"$x / \lambda_{D}$", fontsize=20)
        ax.set_ylabel(r"$y / \lambda_{D}$", fontsize=20)
        ax.set_ylim(0, d["gp"]-1)
        ax.set_xlim(0, d["gp"]-1)
        ax.set_aspect("equal")
        pyplot.savefig(pp1, format="pdf")
        pyplot.close()

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", 0.1*j*(d["steps"]/n)), fontsize=25, y=0.99)
        colors = ax.pcolormesh(x, y, Ey, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$E_{y} \left(x,y \right) \left(e / m_{e} \omega_p v_{th} \right)$", fontsize=20)
        colors.set_clim(numpy.min(Ey_only), numpy.max(Ey_only))
        ax.set_xlabel(r"$x / \lambda_{D}$", fontsize=20)
        ax.set_ylabel(r"$y / \lambda_{D}$", fontsize=20)
        ax.set_ylim(0, d["gp"]-1)
        ax.set_xlim(0, d["gp"]-1)
        ax.set_aspect("equal")
        pyplot.savefig(pp2, format="pdf")
        pyplot.close()

    pp1.close()
    pp2.close()
    pp.close()

def plot_phi(file, n):
    with open(file) as json_data:
        d = json.load(json_data)
    
    Potential = []
    to_plot = [step for step in numpy.arange(0, d["steps"], d["steps"] / n, dtype=int)]  

    for t in to_plot:
        data = numpy.loadtxt("{}/phi/step{}.dat".format(d["Directory"], t), unpack=True)
        Potential.append(data)

    pp = PdfPages("{}/Electric_Potential.pdf".format(d["Directory"]))
    phi_only = []

    for j in range(len(Potential)):
        phi_only.append(Potential[j][2])
    phi_only = numpy.array(phi_only)

    for j in tqdm(range(len(Potential))):
        x = Potential[j][0].reshape((d["gp"], d["gp"]))
        y = Potential[j][1].reshape((d["gp"], d["gp"]))
        phi = Potential[j][2].reshape((d["gp"], d["gp"]))

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", 0.1*j*(d["steps"]/n)), fontsize=25, y=0.99)
        colors = ax.pcolormesh(x, y, phi, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$\phi \left(x,y \right) \left(e / kT \right)$", fontsize=20)
        colors.set_clim(numpy.min(phi_only), numpy.max(phi_only))
        ax.set_xlabel(r"$x / \lambda_{D}$", fontsize=20)
        ax.set_ylabel(r"$y / \lambda_{D}$", fontsize=20)
        ax.set_ylim(0, d["gp"]-1)
        ax.set_xlim(0, d["gp"]-1)
        ax.set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

    pp.close()   

def plot_rho(file, n):
    with open(file) as json_data:
        d = json.load(json_data)
    density = []

    to_plot = [step for step in numpy.arange(0, d["steps"], d["steps"]/n, dtype=int)]  

    for t in to_plot:
        data = numpy.loadtxt("{}/rho/step{}.dat".format(d["Directory"], t), unpack=True)
        density.append(data)

    pp = PdfPages("{}/chargeDensity.pdf".format(d["Directory"]))
    rho_only = []

    for j in range(len(density)):
        rho_only.append(density[j][2])
    rho_only = numpy.array(rho_only)

    for j in tqdm(range(len(density))):
        x = density[j][0].reshape((d["gp"], d["gp"]))
        y = density[j][1].reshape((d["gp"], d["gp"]))
        rho = density[j][2].reshape((d["gp"], d["gp"]))

        fig = pyplot.figure()
        ax = fig.add_subplot(111)
        pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", 0.1*j*(d["steps"]/n)), fontsize=25, y=0.99)
        colors = ax.pcolormesh(x, y, rho, cmap="jet", shading="gouraud")
        bar = fig.colorbar(colors, ax=ax)
        bar.set_label(r"$\rho(x, y)$", fontsize=20)
        colors.set_clim(numpy.min(rho_only), numpy.max(rho_only))
        ax.set_xlabel(r"$x / \lambda_{D}$", fontsize=20)
        ax.set_ylabel(r"$y / \lambda_{D}$", fontsize=20)
        ax.set_ylim(0, d["gp"]-1)
        ax.set_xlim(0, d["gp"]-1)
        ax.set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()
        
    pp.close()    

@click.command()
@click.option("--n", type=int)
@click.argument("file")
def main(file, n):
    plot_Efield(file, n)
    plot_phi(file, n)
    plot_rho(file, n)

if __name__ == '__main__':
    main()