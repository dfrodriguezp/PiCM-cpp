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

    to_plot = [i for i in range(0, d["steps"], 200)]  
    
    vels_x = []
    vels_y = []

    for i in to_plot:
        data = numpy.loadtxt(d["directory"] + "/velocities/step" + str(i) + ".dat", unpack=True)   
        vels_x.append(data[0])
        vels_y.append(data[1])

    pp = PdfPages(d["directory"] + "/velocities.pdf")

    for j in range(len(vels_x)):
        pyplot.figure()
        pyplot.hist(vels_x[j], bins=20, color="forestgreen", label="step " + str(to_plot[j]), normed=True)
        pyplot.xlabel("$v_x$ [au]", size=20)
        pyplot.ylabel("$f(v_x)$", size=20)
        pyplot.legend(loc=0)
        pyplot.tight_layout()
        pyplot.grid()
        pyplot.savefig(pp, format="pdf")

        pyplot.figure()
        pyplot.hist(vels_y[j], bins=20, color="forestgreen", label="step " + str(to_plot[j]), normed=True)
        pyplot.xlabel("$v_y$ [au]", size=20)
        pyplot.ylabel("$f(v_y)$", size=20)
        pyplot.legend(loc=0)
        pyplot.tight_layout()
        pyplot.grid()
        pyplot.savefig(pp, format="pdf")        
        pyplot.close()

    pp.close()
if __name__ == '__main__':
    plot()