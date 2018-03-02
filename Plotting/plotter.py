import numpy
import json
from matplotlib import pyplot
from matplotlib import animation
import click

@click.command()
@click.option("--name", type=str)
@click.option("--x", type=str)
@click.option("--y", type=str)
@click.argument("file")
def plotter(file, name, x, y):
    with open(file) as json_data:
        d = json.load(json_data)

    fig = pyplot.figure(figsize=(16,9))
    ax = pyplot.axes(xlim=(0.0, d["box size"]), ylim=(-0.5, 0.5))
    scat = ax.scatter([], [], color="0.2", lw=0, s=10)
    ax.set_xlabel("${}$ [au]".format(x), size=20)
    ax.set_ylabel("${}$ [au]".format(y), size=20)

    def init():
        scat.set_offsets([])
        return scat,

    def animate(i):
        filename = "{}/{}/step{}.dat".format(d["directory"], name, i)
        x, y = numpy.loadtxt(filename, unpack=True)
        data = numpy.hstack((x[:len(x),numpy.newaxis],y[:len(x),numpy.newaxis]))
        scat.set_offsets(data)
        return scat,
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=d["steps"], interval=1, blit=True)
    pyplot.show()
    # anim.save(d["directory"] + "/phase_space_x.mp4", fps=30, dpi=100)

# def plotter(file):
#     with open(file) as json_data:
#         d = json.load(json_data)

#     PS = []
#     to_plot = [1800]

#     for i in to_plot:
#         data = numpy.loadtxt(d['directory'] + '/phase_space_x/step' + str(i) + '.dat', unpack=True)
#         PS.append(data)

#     for j in range(len(PS)):
#         pyplot.figure()
#         pyplot.scatter(PS[j][0], PS[j][1], color='k', lw=0, s=10)
#         pyplot.xlabel('$x(t)$ [au]', size=20)
#         pyplot.ylabel('$v(t)$ [au]', size=20)
#         pyplot.xlim(0, d['size'])
#         pyplot.ylim(-0.6, 0.6)
#         # pyplot.ylim(min(PS[j][1]), max(PS[j][1]))
#         pyplot.savefig(d['directory'] + '/phase_space_x_' + str(to_plot[j]) + '.png')
#         pyplot.close()

if __name__ == '__main__':
    plotter()
