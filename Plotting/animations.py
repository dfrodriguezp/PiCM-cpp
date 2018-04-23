import numpy
import json
from matplotlib import pyplot
from matplotlib import animation
from matplotlib.backends.backend_pdf import PdfPages
import click

@click.command()
@click.option("--name", type=str)
@click.argument("file")
def main(file, name):
    with open(file) as json_data:
        d = json.load(json_data)

    fig = pyplot.figure(figsize=(16,9))
    if name == "phaseSpace":
        ax = pyplot.axes(xlim=(0.0, d["Box size"]), ylim=(-2.5, 2.5))
        ax.set_xlabel("$x$ [au]", size=20)
        ax.set_ylabel("$v_x$ [au]", size=20)
    elif name == "space":
        ax = pyplot.axes(xlim=(0.0, d["Box size"]), ylim=(0.0, d["Box size"]))
        ax.set_xlabel("$x$ [au]", size=20)
        ax.set_ylabel("$y$ [au]", size=20)
        ax.set_aspect("equal")
    elif name == "velocities":
        ax = pyplot.axes(xlim=(-0.1, 0.1), ylim=(-1.0, 1.0))
        ax.set_xlabel("$v_x$ [au]", size=20)
        ax.set_ylabel("$v_y$ [au]", size=20) 

    scat = ax.scatter([], [], color="black", lw=0, s=10)

    def init():
        scat.set_offsets([])
        return scat,

    def animate(i):
        filename = "{}/{}/step{}.dat".format(d["Directory"], name, i)
        x, y = numpy.loadtxt(filename, unpack=True)
        data = numpy.hstack((x[:len(x),numpy.newaxis],y[:len(x),numpy.newaxis]))
        scat.set_offsets(data)
        return scat,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=d["Steps"], interval=1, blit=True)
    pyplot.show()
    # anim.save(d["directory"] + "/phase_space_x.mp4", fps=30, dpi=100)

if __name__ == '__main__':
    main()