from matplotlib import pyplot
import numpy
import json
import os
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
        print("ERROR! You want to plot a snapshot {} but the limit is {}.".format(
            step, root["steps"]-1))
        return

    output = root["output"]

    # Initial data for colors
    sign = numpy.loadtxt(root["sample"], usecols=(5,))

    colors = numpy.full(len(sign), "r")
    colors[sign < 0] = "b"

    if "phase_space" in root["results"]:
        try:
            x, y, _, _, _ = numpy.loadtxt(
                "{}/phase_space/step_{}_.dat".format(output, step), unpack=True)
        except OSError as error:
            print("ERROR! You want to plot a snapshot that doesn't exist.")
            return
    else:
        try:
            x, y = numpy.loadtxt(
                "{}/space/step_{}_.dat".format(output, step), unpack=True)
        except OSError as error:
            print("ERROR! You want to plot a snapshot that doesn't exist.")
            return

    if len(x) > 1e5 or numpy.allclose(len(x), 1e5, atol=500):
        colors = colors[::100]
        x = x[::100]
        y = y[::100]

    elif (1e5 > len(x) > 0.5e5) or numpy.allclose(len(x), 0.5e5, atol=500):
        colors = colors[::50]
        x = x[::50]
        y = y[::50]

    Lx, Ly = root["sys_length"]
    t = root["dt"] * step

    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.scatter(x, y, c=colors)
    pyplot.xlabel(r"$x / \lambda_D$", fontsize=25)
    pyplot.ylabel(r"$y / \lambda_D$", fontsize=25)
    pyplot.xlim(0, Lx)
    pyplot.ylim(0, Ly)
    if Lx == Ly:
        pyplot.gca().set_aspect("equal")
    pyplot.grid()
    if not "space" in root["results"]:
        os.system("mkdir -p {}/space".format(output))
    pyplot.savefig("{}/space/step_{}_.pdf".format(output, step))
    pyplot.close()


if __name__ == "__main__":
    main()
