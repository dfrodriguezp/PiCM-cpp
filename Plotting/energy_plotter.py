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
    pyplot.plot(t*d["dt"], KE, "--o", ms=2)
    pyplot.axvline(x=t[numpy.argmin(KE)]*d["dt"], ls="--", color="red")
    pyplot.ylabel(r"$E_{k} \left(\varepsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    pyplot.xticks(size=20)
    pyplot.yticks(size=20)
    pyplot.grid()
    pyplot.savefig("{}/Kinetic Energy.pdf".format(d["Directory"]))
    pyplot.close()

    pyplot.figure(figsize=(12, 9))
    pyplot.suptitle("Electric Field Energy", fontsize=25)
    pyplot.plot(t*d["dt"], FE, "--o", ms=2)
    pyplot.axvline(x=t[numpy.argmax(FE)]*d["dt"], ls="--", color="red")
    pyplot.ylabel(r"$E_F \left(\varepsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    pyplot.xticks(size=20)
    pyplot.yticks(size=20)
    pyplot.grid()
    pyplot.savefig("{}/Field Energy.pdf".format(d["Directory"]))
    pyplot.close()

if __name__ == '__main__':
    main()