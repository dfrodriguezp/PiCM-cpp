import numpy
from matplotlib import pyplot
import click
import json
import h5py

@click.command()
@click.argument("file", default="../TSI_test.h5")
def main(file):
    f = h5py.File(file, "r")
    energy_data = f["energy"]

    t, KE, FE = range(len(energy_data)), energy_data[:,0], energy_data[:,1]

    pyplot.figure(figsize=(12, 9))
    pyplot.suptitle("Kinetic Energy", fontsize=25)
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
    pyplot.suptitle("Electric Field Energy", fontsize=25)
    pyplot.plot(t, FE, "--o", ms=2)
    # pyplot.axvline(x=t[numpy.argmax(FE)]*d["dt"], ls="--", color="red")
    pyplot.ylabel(r"$E_F \left(\varepsilon_{0} / n_{0}kT \right)$", fontsize=35)
    pyplot.xlabel(r"$\omega_{pe}t$", fontsize=35)
    pyplot.xticks(size=20)
    pyplot.yticks(size=20)
    pyplot.grid()
    pyplot.savefig("Field_Energy.pdf")
    pyplot.close()

if __name__ == '__main__':
    main()