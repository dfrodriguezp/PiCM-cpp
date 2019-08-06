from matplotlib import pyplot
import numpy
import json
import os
import click
import matplotlib
matplotlib.style.use('classic')


@click.command()
@click.argument("jsonfile")
@click.argument("step", default=0)
def main(jsonfile, step):
    with open(jsonfile, "r") as json_input:
        root = json.load(json_input)

    if step >= root["steps"]:
        print("ERROR! You want to plot a snapshot {} but the limit is {}.".format(
            step, root["steps"]-1))
        return

    output = root["output"]
    if "phase_space" in root["results"]:
        try:
            _, _, vx, vy, vz = numpy.loadtxt(
                "{}/phase_space/step_{}_.dat".format(output, step), unpack=True)
        except OSError as error:
            print("ERROR! You want to plot a snapshot that doesn't exist.")
            return
    else:
        try:
            vx, vy, vz = numpy.loadtxt(
                "{}/velocities/step_{}_.dat".format(output, step), unpack=True)
        except OSError as error:
            print("ERROR! You want to plot a snapshot that doesn't exist.")
            return

    if not "velocities" in root["results"]:
        os.system("mkdir -p {}/velocities".format(output))

    vtotal = numpy.linalg.norm([vx, vy, vz], axis=0)
    t = root["dt"] * step
    # x-component
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.hist(vx, bins=50, color="red")
    pyplot.xlabel(r"$v_x / v_{\rm{th}}$", fontsize=25)
    pyplot.grid()
    pyplot.savefig("{}/velocities/step_{}_vx_.pdf".format(output, step))
    pyplot.close()

    # y-component
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.hist(vy, bins=50, color="green")
    pyplot.xlabel(r"$v_y / v_{\rm{th}}$", fontsize=25)
    pyplot.grid()
    pyplot.savefig("{}/velocities/step_{}_vy_.pdf".format(output, step))
    pyplot.close()

    # z-component
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.hist(vz, bins=50, color="blue")
    pyplot.xlabel(r"$v_z / v_{\rm{th}}$", fontsize=25)
    pyplot.grid()
    pyplot.savefig("{}/velocities/step_{}_vz_.pdf".format(output, step))
    pyplot.close()

    # v total
    pyplot.figure()
    pyplot.title(r"$\omega_{\rm{pe}}$" + f"$t = {t}$", fontsize=25)
    pyplot.hist(vtotal, bins=50, color="yellow")
    pyplot.xlabel(r"$\left|\vec{v}\right| / v_{\rm{th}}$", fontsize=25)
    pyplot.grid()
    pyplot.savefig("{}/velocities/step_{}_vtotal_.pdf".format(output, step))
    pyplot.close()


if __name__ == "__main__":
    main()
