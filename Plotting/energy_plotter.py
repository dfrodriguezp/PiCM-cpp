import numpy
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import click
import json

@click.command()
@click.argument("file")
def plot(file):
    with open(file) as json_data:
        d = json.load(json_data)
    Ek, EF, t = numpy.loadtxt(d["directory"] + "/energies/energies.dat", unpack=True)
    
    # pp = PdfPages(d["directory"] + "/energies.pdf")

    pyplot.figure()
    pyplot.plot(t[:1000], Ek[:1000]/max(Ek[:1000]), "--o", ms=2)
    pyplot.ylabel("$E_k$ [au]", size=20)
    pyplot.xlabel("$Steps$", size=20)
    pyplot.tight_layout()
    pyplot.grid()
    pyplot.savefig("%s/E_k.jpg" %(d["directory"]))
    # pyplot.savefig(pp, format="pdf")

    pyplot.figure()
    pyplot.plot(t[:1000], EF[:1000]/max(EF[:1000]), "--o", ms=2)
    pyplot.ylabel("$E_F$ [au]", size=20)
    pyplot.xlabel("$Steps$", size=20)
    pyplot.tight_layout()
    pyplot.grid()
    pyplot.savefig("%s/E_F.jpg" %(d["directory"]))
    # pyplot.savefig(pp, format="pdf")

    pyplot.close()
    # pp.close()
if __name__ == '__main__':
    plot()