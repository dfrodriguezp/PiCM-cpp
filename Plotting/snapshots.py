import numpy
import json
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import click
from tqdm import tqdm

def downloadData(name, toPlot, jsData):
    data = list()
    for i in toPlot:
        download = numpy.loadtxt("{}/{}/step{}.dat".format(jsData["Directory"], name, i), unpack=True)
        data.append(download)
    pp = PdfPages("{}/{}.pdf".format(jsData["Directory"], name))

    return data, pp

@click.command()
@click.option("--name", type=str)
@click.option("--n", type=int)
@click.argument("file")
def main(file, name, n):
    with open(file) as json_data:
        d = json.load(json_data)

    toPlot = [step for step in numpy.arange(0, d["steps"], d["steps"] / n, dtype=int)]

    names = ["phaseSpace", "space", "velocities"]

    if name in names:
        data, pp = downloadData(name, toPlot, d)

    elif name == "directions":
        pos = list()
        vel, pp = downloadData("velocities", toPlot, d)
        pos, pp = downloadData("space", toPlot, d)
        pp = PdfPages("{}/directions.pdf".format(d["Directory"]))       

    for j in tqdm(range(len(data))):
        pyplot.figure()
        pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", 0.1*j*(d["steps"]/n)), fontsize=30, y=0.99)
        if name == "phaseSpace":  
            pyplot.scatter(data[j][0][::100], data[j][1][::100], color="black", lw=0, s=5)
            pyplot.xlabel(r"$x / \lambda_{D}$", fontsize=20)
            pyplot.ylabel(r"$v_x / v_{th}$", fontsize=20)
            pyplot.xlim(0, d["L"])
            pyplot.ylim(-15, 15)
        elif name == "space":
            pyplot.scatter(data[j][0], data[j][1], color="black", lw=0, s=5)      
            pyplot.xlabel(r"$x / \lambda_{D}$", fontsize=20)
            pyplot.ylabel(r"$y / \lambda_{D}$", fontsize=20)
            pyplot.xlim(0, d["L"])
            pyplot.ylim(0, d["L"])
            pyplot.gca().set_aspect("equal")
        elif name == "velocities":
            pyplot.scatter(data[j][0], data[j][1], color="black", lw=0, s=5)      
            pyplot.xlabel(r"$v_x / v_{th}$", fontsize=20)
            pyplot.ylabel(r"$v_y / v_{th}$", fontsize=20)
            pyplot.xlim(-0.1, 0.1)
            pyplot.ylim(-0.1, 0.1)
        elif name == "directions":
            pyplot.quiver(pos[j][0], pos[j][1], vel[j][0], vel[j][0], pivot="tail")      
            pyplot.xlabel(r"$v_x / v_{th}$", fontsize=20)
            pyplot.ylabel(r"$v_y / v_{th}$", fontsize=20)
            pyplot.xlim(0, d["L"])
            pyplot.ylim(0, d["L"])
            pyplot.gca().set_aspect("equal")
        pyplot.savefig(pp, format="pdf")
        pyplot.close()

    pp.close()

if __name__ == '__main__':
    main()