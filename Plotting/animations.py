import numpy
import json
from matplotlib import pyplot
from matplotlib import animation
from matplotlib.backends.backend_pdf import PdfPages
import click

@click.command()
@click.option("--name", type=str)
@click.option("--vlimx", type=float)
@click.option("--vlimy", type=float)
@click.argument("file")
def main(file, name):
    with open(file) as json_data:
        d = json.load(json_data)

    fig = pyplot.figure(figsize=(16,9))
    if name == "phaseSpace":
        ax = pyplot.axes(xlim=(0.0, d["L"]), ylim=(-vlimx, vlimx))
        ax.set_xlabel("$x$ [au]", size=20)
        ax.set_ylabel("$v_x$ [au]", size=20)
    elif name == "space":
        ax = pyplot.axes(xlim=(0.0, d["L"]), ylim=(0.0, d["L"]))
        ax.set_xlabel("$x$ [au]", size=20)
        ax.set_ylabel("$y$ [au]", size=20)
        ax.set_aspect("equal")
    elif name == "velocities":
        ax = pyplot.axes(xlim=(-vlimx, vlimx), ylim=(-vlimy, vlimy))
        ax.set_xlabel("$v_x$ [au]", size=20)
        ax.set_ylabel("$v_y$ [au]", size=20) 

    scat = ax.scatter([], [], color="black", lw=0, s=10)

    def init():
        scat.set_offsets([])
        return scat,

    def animate(i):
        filename = "{}/{}/step{}.dat".format(d["Directory"], name, i)
        x, y = numpy.loadtxt(filename, unpack=True)
        data = numpy.hstack((x[:len(x),numpy.newaxis][::200],y[:len(x),numpy.newaxis][::200]))
        scat.set_offsets(data)
        return scat,

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=d["steps"], interval=1, blit=True)
    pyplot.show()
    # print("Saving...")
    # anim.save(d["Directory"] + "/phase_space_x.mp4", fps=30, dpi=100)
    # print("Saved.")

if __name__ == '__main__':
    main()