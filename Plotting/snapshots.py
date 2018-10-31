import numpy
import json
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages
import click
from tqdm import tqdm

def downloadData(name, toPlot, jsData):
    data = list()
    colors = list()
    for i in toPlot:
        download = numpy.loadtxt("{}/{}/step{}.dat".format(jsData["Directory"], name, i), unpack=True)
        data.append(download)
    pp = PdfPages("{}/{}.pdf".format(jsData["Directory"], name))
    # pp = "something"
    N = len(data[0][1])
    for i in range(N):
        if i < (N/2):
            colors.append("blue")
        else:
            colors.append("red")
    return data, pp, colors

@click.command()
@click.option("--name", type=str)
@click.option("--n", type=int)
@click.option("--vlimx", type=float)
@click.option("--vlimy", type=float)
@click.option("--plot", default=True)
@click.argument("file")
def main(file, name, n, vlimx, vlimy, plot):
    with open(file) as json_data:
        d = json.load(json_data)

    toPlot = [step for step in numpy.arange(0, d["steps"], d["steps"] / n, dtype=int)]
    # toPlot = numpy.array([0, 13, 22.5, 27, 35, 80]) / d["dt"]
    # toPlot = numpy.array(toPlot, dtype=int)
    names = ["phaseSpace", "space", "velocities"]

    if name in names:
        data, pp, colors = downloadData(name, toPlot, d)
        N = len(data)
    elif name == "directions":
        pos = list()
        vel, pp, colors = downloadData("velocities", toPlot, d)
        pos, pp, colors = downloadData("space", toPlot, d)
        pp = PdfPages("{}/directions.pdf".format(d["Directory"]))       
        N = len(pos)
    if plot:
        for j in tqdm(range(N)):
            pyplot.figure()
            pyplot.suptitle(r"$\omega_{}t={:.1f}$".format("{pe}", d["dt"]*toPlot[j]), fontsize=30, y=0.99)
            if name == "phaseSpace":  
                pyplot.scatter(data[j][0][::100], data[j][1][::100], color=colors[::100], lw=0, s=5)
                pyplot.xlabel(r"$x / \lambda_{D}$", fontsize=20)
                pyplot.ylabel(r"$v_x / v_{th}$", fontsize=20)
                pyplot.xlim(0, d["L"])
                pyplot.ylim(-vlimx, vlimx)
            elif name == "space":
                pyplot.scatter(data[j][0][::200], data[j][1][::200], color=colors[::200], lw=0, s=5)      
                pyplot.xlabel(r"$x / \lambda_{D}$", fontsize=20)
                pyplot.ylabel(r"$y / \lambda_{D}$", fontsize=20)
                pyplot.xlim(0, d["L"])
                pyplot.ylim(0, d["L"])
                pyplot.gca().set_aspect("equal")
            elif name == "velocities":
                pyplot.scatter(data[j][0][::200], data[j][1][::200], color=colors[::200], lw=0, s=5)      
                pyplot.xlabel(r"$v_x / v_{th}$", fontsize=20)
                pyplot.ylabel(r"$v_y / v_{th}$", fontsize=20)
                pyplot.xlim(-vlimx, vlimx)
                pyplot.ylim(-vlimy, vlimy)
            elif name == "directions":
                pyplot.quiver(pos[j][0][::1000], pos[j][1][::1000], vel[j][0][::1000], vel[j][0][::1000], color=colors[::1000], pivot="tail")      
                pyplot.xlabel(r"$x / \lambda_{D}$", fontsize=20)
                pyplot.ylabel(r"$y / \lambda_{D}$", fontsize=20)
                pyplot.xlim(0, d["L"])
                pyplot.ylim(0, d["L"])
                pyplot.gca().set_aspect("equal")
            # pyplot.savefig("{}/{}_wpt_{}.png".format(d["Directory"], name, d["dt"]*toPlot[j]), dpi=300)
            pyplot.savefig(pp, format="pdf")
            pyplot.close()

        pp.close()

if __name__ == '__main__':
    main()