import numpy
import json
import click

@click.command()
@click.option("-gp", default=8, type=int)
@click.option("-system", default="two_stream", type=str)
@click.option("-np", default=200, type=int)
@click.option("-dr", default=1.0, type=float)
@click.option("-dt", default=0.1, type=float)
@click.option("-vt", default=1.0, type=float)
@click.option("-vd", default=2.0, type=float)
@click.option("-extb", default=(0.0, 0.0, 0.0))
@click.option("-steps", default=100, type=int)
@click.option("-seed", default=numpy.random.randint(10000, 1000000))
def main(gp, system, np, dr, dt, vt, vd, extb, steps, seed):
    numpy.random.seed(seed)
    Nm = int(np / 2)
    L = (gp - 1) * dr
    electronic_density = np / (L*L)

    if system == "two_stream":
        pos, vel = two_stream(L, gp, Nm, dr, vt, vd)

        datafile = "system_{}.dat".format(system)

        output = open(datafile, "w")
        for i in range(int(len(vel)/2)):
            output.write("{:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} -1.0 1\n".format(*pos[i], *vel[i], electronic_density))
        for i in range(int(len(vel)/2), len(vel)):
            output.write("{:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} 1.0 0\n".format(*pos[i], *vel[i], electronic_density))
        output.close()

        jsondata = {
            "steps": steps,
            "seed": seed,
            "dr": dr,
            "npoints": gp,
            "dt": dt,
            "N": np,
            "vt": vt,
            "vd": vd,
            "B": {
                "Bx": extb[0],
                "By": extb[1],
                "Bz": extb[2]
                },
            "output": datafile
        }

        with open("parameters.json", "w") as file:
            json.dump(jsondata, file)

def two_stream(L, gp, N, dr, vt, vd):
    pos = list()
    vel = list()
    spacing = numpy.linspace(dr, L-dr, int(numpy.sqrt(N)))
    for x in spacing:
        for y in spacing:
            pos.append([x, y, 0.0])

    indexes = numpy.random.choice(range(N), N, replace=False)
    right = indexes[:int(N/2)]
    left = indexes[int(N/2):]

    pos_write = list()
    for i in right:
        vel_right = numpy.random.normal(vd, vt)
        vel.append([vel_right, 0.0, 0.0])
        pos_write.append(pos[i])
    for i in left:
        vel_left = numpy.random.normal(-vd, vt)
        vel.append([vel_left, 0.0, 0.0])
        pos_write.append(pos[i])

    for x in spacing:
        for y in spacing:
            pos_write.append([x, y, 0.0])

    for i in range(N):
        vel.append([0.0, 0.0, 0.0])

    pos = numpy.array(pos_write)
    vel = numpy.array(vel)

    return pos, vel 

if __name__ == '__main__':
    main()