import numpy
from matplotlib import pyplot
import click
import json

@click.command()
@click.argument("file", default="../results/energy/energy.dat")
def main(file):
    t, KE, FE = numpy.loadtxt(file, unpack=True)

    pyplot.figure(figsize=(12, 9))
    pyplot.plot(t, KE, "--o", ms=2)
    # pyplot.axvline(x=t[numpy.argmin(KE)]*d["dt"], ls="--", color="red")
    pyplot.ylabel(r"$E_{k} \left(\varepsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    pyplot.xticks(size=20)
    pyplot.yticks(size=20)
    pyplot.grid()
    pyplot.savefig("Kinetic_Energy.pdf")
    pyplot.close()

    pyplot.figure(figsize=(12, 9))
    pyplot.plot(t, FE, "--o", ms=2)
    # pyplot.axvline(x=t[numpy.argmax(FE)]*d["dt"], ls="--", color="red")
    pyplot.ylabel(r"$E_F \left(\varepsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    # pyplot.ylim(0, 1000)
    pyplot.xticks(size=20)
    pyplot.yticks(size=20)
    pyplot.grid()
    pyplot.savefig("Field_Energy.pdf")
    pyplot.close()

if __name__ == '__main__':
    main()