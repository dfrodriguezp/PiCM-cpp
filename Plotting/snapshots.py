import numpy
import json
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import click
from tqdm import tqdm

@click.command()
@click.option("--name", type=str)
@click.option("--n", type=int)
@click.argument("file")
def main(file, name, n):
    with open(file) as json_data:
        d = json.load(json_data)

    data = list()
    to_plot = [step for step in numpy.arange(0, d["Steps"], d["Steps"] / n, dtype=int)]

    if name == "phaseSpace":
        for i in to_plot:
            download = numpy.loadtxt("{}/phaseSpace/step{}.dat".format(d["Directory"], i), unpack=True)
            data.append(download)
        pp = PdfPages("{}/phaseSpace.pdf".format(d["Directory"]))

    elif name == "space":
        for i in to_plot:
            download = numpy.loadtxt("{}/space/step{}.dat".format(d["Directory"], i), unpack=True)
            data.append(download)
        pp = PdfPages("{}/space.pdf".format(d["Directory"]))

    elif name == "velocities":
        for i in to_plot:
            download = numpy.loadtxt("{}/velocities/step{}.dat".format(d["Directory"], i), unpack=True)
            data.append(download)
        pp = PdfPages("{}/velocities.pdf".format(d["Directory"]))

    elif name == "directions":
        pos = list()
        for i in to_plot:
            download = numpy.loadtxt("{}/velocities/step{}.dat".format(d["Directory"], i), unpack=True)
            download2 = numpy.loadtxt("{}/space/step{}.dat".format(d["Directory"], i), unpack=True)
            data.append(download)
            pos.append(download2)

        pp = PdfPages("{}/directions.pdf".format(d["Directory"]))       

    for j in tqdm(range(len(data))):
        pyplot.figure()
        pyplot.suptitle(r"$t={:.1f} \ \omega_p^{}$".format(0.1*j*(d["Steps"]/n), "{-1}"), fontsize=30, y=0.99)
        if name == "phaseSpace":  
            pyplot.scatter(data[j][0], data[j][1], color="black", lw=0, s=5)
            pyplot.xlabel("$x$ [au]", size=20)
            pyplot.ylabel("$v_x$ [au]", size=20)
            pyplot.xlim(0, d["Box size"])
            pyplot.ylim(-0.1, 0.1)
        elif name == "space":
            pyplot.scatter(data[j][0], data[j][1], color="black", lw=0, s=5)      
            pyplot.xlabel("$x$ [au]", size=20)
            pyplot.ylabel("$y$ [au]", size=20)
            pyplot.xlim(0, d["Box size"])
            pyplot.ylim(0, d["Box size"])
            pyplot.gca().set_aspect("equal")
        elif name == "velocities":
            pyplot.scatter(data[j][0], data[j][1], color="black", lw=0, s=5)      
            pyplot.xlabel("$v_x$ [au]", size=20)
            pyplot.ylabel("$v_y$ [au]", size=20)
            pyplot.xlim(-0.1, 0.1)
            pyplot.ylim(-0.1, 0.1)
        elif name == "directions":
            pyplot.quiver(pos[j][0], pos[j][1], data[j][0], data[j][0], pivot="tail")      
            pyplot.xlabel("$v_x$ [au]", size=20)
            pyplot.ylabel("$v_y$ [au]", size=20)
            pyplot.xlim(0, d["Box size"])
            pyplot.ylim(0, d["Box size"])
            pyplot.gca().set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

    pp.close()

if __name__ == '__main__':
    main()