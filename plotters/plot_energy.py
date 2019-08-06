import numpy
import json
import click
from matplotlib import pyplot
from matplotlib import style
style.use("classic")


@click.command()
@click.argument("jsonfile")
def main(jsonfile):
    with open(jsonfile, "r") as json_input:
        root = json.load(json_input)

    output = root["output"]
    dt = root["dt"]
    step, KE, FE = numpy.loadtxt(
        "{}/energy/energy.dat".format(output), unpack=True)

    total = KE + FE
    t = step * dt

    pyplot.figure()
    pyplot.plot(t, KE, label="Kinetic")
    pyplot.plot(t, FE, label="Field")
    pyplot.plot(t, total, label="Total")
    pyplot.xlabel(r"$\omega_{\rm{pe}}t$", fontsize=25)
    pyplot.ylabel(r"$E / (n_0 T_e / \varepsilon_0)$", fontsize=25)
    pyplot.grid()
    pyplot.legend(loc="best")
    pyplot.tight_layout()
    pyplot.savefig("{}/energy/energy.pdf".format(output))
    pyplot.close()


if __name__ == "__main__":
    main()
