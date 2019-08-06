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
        x, y, Ex, Ey = numpy.loadtxt(
            "{}/Efield/step_{}_.dat".format(output, step), unpack=True)
    except OSError as error:
        print("ERROR! You want to plot a snapshot that doesn't exist.")
        return

    Etotal = numpy.linalg.norm([Ex, Ey], axis=0)
    Nx, Ny = root["grid_size"]
    Lx, Ly = root["sys_length"]
    dx, dy = Lx / Nx, Ly / Ny
    x = x.reshape((Nx, Ny))
    y = y.reshape((Nx, Ny))
    Ex = Ex.reshape((Nx, Ny))
    Ey = Ey.reshape((Nx, Ny))
    Etotal = Etotal.reshape((Nx, Ny))
    t = root["dt"] * step

    # x-component
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    color_map = pyplot.pcolormesh(x, y, Ex, shading="gouraud", cmap="jet")
    bar = pyplot.colorbar(color_map, ax=pyplot.gca())
    pyplot.xlim(0, Lx - dx)
    pyplot.ylim(0, Ly - dy)
    pyplot.xlabel(r"$x / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$y / \lambda_D$", fontsize=25)
    bar.set_label(r"$E_x / (m_e\omega_{\rm{pe}}v_{\rm{th}} / e)$", fontsize=25)
    pyplot.gca().set_aspect("equal")
    pyplot.savefig("{}/Efield/step_{}_Ex_.pdf".format(output, step))
    pyplot.close()

    # y-component
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    color_map = pyplot.pcolormesh(x, y, Ey, shading="gouraud", cmap="jet")
    bar = pyplot.colorbar(color_map, ax=pyplot.gca())
    pyplot.xlim(0, Lx - dx)
    pyplot.ylim(0, Ly - dy)
    pyplot.xlabel(r"$x / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$y / \lambda_D$", fontsize=25)
    bar.set_label(r"$E_y / (m_e\omega_{\rm{pe}}v_{\rm{th}} / e)$", fontsize=25)
    pyplot.gca().set_aspect("equal")
    pyplot.savefig("{}/Efield/step_{}_Ey_.pdf".format(output, step))
    pyplot.close()

    # E total
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    color_map = pyplot.pcolormesh(x, y, Etotal, shading="gouraud", cmap="jet")
    bar = pyplot.colorbar(color_map, ax=pyplot.gca())
    pyplot.xlim(0, Lx - dx)
    pyplot.ylim(0, Ly - dy)
    pyplot.xlabel(r"$x / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$y / \lambda_D$", fontsize=25)
    bar.set_label(
        r"$\left|\vec{E}\right| / (m_e\omega_{\rm{pe}}v_{\rm{th}} / e)$", fontsize=25)
    pyplot.gca().set_aspect("equal")
    pyplot.savefig("{}/Efield/step_{}_Etotal_.pdf".format(output, step))
    pyplot.close()


if __name__ == "__main__":
    main()
