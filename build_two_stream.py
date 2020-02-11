import numpy
from itertools import product
import json
import os

N = int(1e5)  # approximate number of particles
Lx, Ly = 64, 64  # size of the system
Nx, Ny = 64, 64  # number of grid points
dx, dy = Lx / Nx, Ly / Ny  # delta_x and delta_y
Bx, By, Bz = 0.0, 0.0, 0.0  # external magnetic field
vd, vt = 5.0, 1.0  # drift and thermal velocities
steps = 500  # total of steps
dt = 0.1  # time step
output = "electrosctatic"  # output folder

margin = dx / 10  # safety margin so the particles don't initialize in the border

line_positions = numpy.linspace(margin, Lx - margin, int(numpy.sqrt(N)))
positions = numpy.array(list(product(line_positions, repeat=2)))
numpy.random.shuffle(positions)
vel_zero = numpy.zeros(int(N / 2))
vel_left = numpy.random.normal(-vd, vt, size=int(N / 4))
vel_right = numpy.random.normal(vd, vt, size=int(N / 4))
velocities = numpy.concatenate((vel_zero, vel_left, vel_right))
QoverM = numpy.concatenate((numpy.ones(int(N / 2)), -numpy.ones(int(N / 2))))
moves = numpy.concatenate((numpy.zeros(int(N / 2)), numpy.ones(int(N / 2))))

sample_name = "two_stream.dat"

if not(os.path.isdir(output)):
    os.mkdir(output)

path_to_sample = f"{output}/{sample_name}"
output_sample = open(path_to_sample, "w")
for i in range(len(positions)):
    output_sample.write("{:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {}\n".format(
        positions[i][0], positions[i][1], velocities[i], 0.0, 0.0, QoverM[i], int(moves[i])))
output_sample.close()

print("Sample built!\n")

json_data = {
    "N": len(positions),
    "steps": steps,
    "grid_size": [Nx, Ny],
    "sys_length": [Lx, Ly],
    "ss_frequency": 50,
    "dt": dt,
    "sample": path_to_sample,
    "output": output,
    "results": ["phase_space"],
    "Bfield": [Bx, By, Bz]
}

with open("{}/sim_{}".format(output, sample_name.replace(".dat", ".json")), "w") as json_out:
    json.dump(json_data, json_out, indent=2, sort_keys=False)

print(".json file built!\n")
print(json_data)
