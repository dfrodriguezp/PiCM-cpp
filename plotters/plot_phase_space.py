import numpy
import json
import click
from matplotlib import pyplot
import matplotlib
matplotlib.style.use('classic')


@click.command()
@click.argument("jsonfile")
@click.argument("step", default=0)
def main(jsonfile, step):
    with open(jsonfile, "r") as json_input:
        root = json.load(json_input)

    if step >= root["steps"]:
        print("ERROR! You want to plot a snapshot {} but the limit is {}.".format(
            step, root["steps"]-1))
        return

    output = root["output"]

    # Step 0 data for colors
    _, _, vx0, vy0, _ = numpy.loadtxt(
        "{}/phase_space/step_0_.dat".format(output), unpack=True)

    colorsx = numpy.full(len(vx0), "r")
    colorsx[vx0 < 0] = "b"

    colorsy = numpy.full(len(vy0), "r")
    colorsy[vy0 < 0] = "b"

    try:
        x, y, vx, vy, _ = numpy.loadtxt(
            "{}/phase_space/step_{}_.dat".format(output, step), unpack=True)
    except OSError as error:
        print("ERROR! You want to plot a snapshot that doesn't exist.")
        return

    if len(vx) > 1e5 or numpy.allclose(len(vx), 1e5, atol=500):
        colorsx = colorsx[::100]
        colorsy = colorsy[::100]
        x = x[::100]
        y = y[::100]
        vx = vx[::100]
        vy = vy[::100]

    elif (1e5 > len(vx) > 0.5e5) or numpy.allclose(len(vx), 0.5e5, atol=500):
        colorsx = colorsx[::50]
        colorsy = colorsy[::50]
        x = x[::50]
        y = y[::50]
        vx = vx[::50]
        vy = vy[::50]

    Lx, Ly = root["sys_length"]
    t = root["dt"] * step

    # x - dimension
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.scatter(x, vx, c=colorsx, edgecolor="none", alpha=0.9)
    pyplot.xlabel(r"$x / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$v_x / v_{\rm{th}}$", fontsize=25)
    pyplot.xlim(0, Lx)
    pyplot.grid()
    pyplot.tight_layout()
    pyplot.savefig("{}/phase_space/step_{}_x_.pdf".format(output, step))
    pyplot.close()

    # y - dimension
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.scatter(y, vy, c=colorsy, edgecolor="none", alpha=0.9)
    pyplot.xlabel(r"$y / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$v_y \ v_{\rm{th}}$", fontsize=25)
    pyplot.xlim(0, Ly)
    pyplot.grid()
    pyplot.tight_layout()
    pyplot.savefig("{}/phase_space/step_{}_y_.pdf".format(output, step))
    pyplot.close()


if __name__ == "__main__":
    main()
