import numpy
from matplotlib import pyplot
import click
import json

@click.command()
@click.argument("file")
def main(file):
    with open(file) as json_data:
        d = json.load(json_data)
    t, KE, FE = numpy.loadtxt(d["Directory"] + "/energy/energy.dat", unpack=True)

    pyplot.figure(figsize=(12, 9))
    pyplot.suptitle("Kinetic Energy", fontsize=25)
    pyplot.plot(t*0.1, KE, "--o", ms=2)
    pyplot.ylabel(r"$E_{k} \left(\epsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    pyplot.grid()
    pyplot.savefig("{}/Kinetic Energy.pdf".format(d["Directory"]))
    pyplot.close()

    pyplot.figure(figsize=(12, 9))
    pyplot.suptitle("Electric Field Energy", fontsize=25)
    pyplot.plot(t*0.1, FE, "--o", ms=2)
    pyplot.ylabel(r"$E_F \left(\epsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    pyplot.grid()
    pyplot.savefig("{}/Field Energy.pdf".format(d["Directory"]))
    pyplot.close()

if __name__ == '__main__':
    main()