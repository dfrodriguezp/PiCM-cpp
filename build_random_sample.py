import numpy
import json
import os


N = int(1e5)  # approximate number of particles
N_electrons = N // 2
N_ions = N // 2
Lx, Ly = 64, 64  # size of the system
Nx, Ny = 64, 64  # number of grid points
dx, dy = Lx / Nx, Ly / Ny  # delta_x and delta_y
Bx, By, Bz = 0.0, 0.0, 0.1  # external magnetic field
steps = 1000  # total of steps
dt = 0.1  # time step
output = "example_random"  # output folder
sample_name = "sample.dat"

if not(os.path.isdir(output)):
    os.mkdir(output)

path_to_sample = f"{output}/{sample_name}"
output_sample = open(path_to_sample, "w")

for i in range(N_electrons):
    x = numpy.random.uniform(0, Lx)
    y = numpy.random.uniform(0, Ly)
    vx = numpy.random.uniform(-0.5, 0.5)
    vy = numpy.random.uniform(-0.5, 0.5)
    vz = numpy.random.uniform(-0.5, 0.5)
    output_sample.write("{} {} {} {} {} {} {}\n".format(
        x, y, vx, vy, vz, -1.0, 1))

for i in range(N_ions):
    x = numpy.random.uniform(0, Lx)
    y = numpy.random.uniform(0, Ly)
    vx = numpy.random.uniform(-0.5, 0.5)
    vy = numpy.random.uniform(-0.5, 0.5)
    vz = numpy.random.uniform(-0.5, 0.5)
    output_sample.write("{} {} {} {} {} {} {}\n".format(
        x, y, vx, vy, vz, 1.0, 1))

output_sample.close()

print("Sample built!\n")

json_data = {
    "N": N,
    "steps": steps,
    "grid_size": [Nx, Ny],
    "sys_length": [Lx, Ly],
    "dt": dt,
    "output": output,
    "sample": path_to_sample,
    "results": ["electric_potential", "space"],
    "Bfield": [Bx, By, Bz]
}

with open("{}/sim_{}".format(output, sample_name.replace(".dat", ".json")), "w") as json_out:
    json.dump(json_data, json_out)

print(".json file built!\n")
print(json_data)
