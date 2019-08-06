from matplotlib import pyplot
import numpy
import json
import click
import matplotlib
matplotlib.style.use('classic')


@click.command()
@click.argument("jsonfile")
@click.argument("step", default=0)
def main(jsonfile, step):
    with open(jsonfile, "r") as json_input:
        root = json.load(json_input)

    if step >= root["steps"]:
        print("ERROR! You want to plot snapshot {} but the limit is {}.".format(
            step, root["steps"]-1))
        return

    output = root["output"]
    try:
        x, y, phi = numpy.loadtxt(
            "{}/phi/step_{}_.dat".format(output, step), unpack=True)
    except OSError as error:
        print("ERROR! You want to plot a snapshot that doesn't exist.")
        return

    Nx, Ny = root["grid_size"]
    Lx, Ly = root["sys_length"]
    dx, dy = Lx / Nx, Ly / Ny
    x = x.reshape((Nx, Ny))
    y = y.reshape((Nx, Ny))
    phi = phi.reshape((Nx, Ny))
    t = root["dt"] * step

    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    color_map = pyplot.pcolormesh(x, y, phi, shading="gouraud", cmap="jet")
    bar = pyplot.colorbar(color_map, ax=pyplot.gca())
    pyplot.xlim(0, Lx - dx)
    pyplot.ylim(0, Ly - dy)
    pyplot.xlabel(r"$x / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$y / \lambda_D$", fontsize=25)
    bar.set_label(r"$\phi / (T_e / e)$", fontsize=25)
    pyplot.gca().set_aspect("equal")
    pyplot.savefig("{}/phi/step_{}_.pdf".format(output, step))
    pyplot.close()


if __name__ == "__main__":
    main()
