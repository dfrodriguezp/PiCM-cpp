import numpy
import json
import glob

def main():
    folder = "random_particles_Bz"
    data = {
    "steps": 2000,
    "seed": 561765,
    "dr": 1.0,
    "gp": 64,
    "dt": 0.05,
    "L": 63.0,
    "N": 504100,
    "vt": 1.0,
    "vd": 5.0,
    "field": {
        "Bx": 0.0,
        "By": 0.0,
        "Bz": 0.1
        }
    }
    outfile = "../{}.json".format(folder)
    data["Directory"] = "/home/daniel/Desktop/random_particles/random_particles_Bz"
    with open(outfile, "w") as file:
        json.dump(data, file)

if __name__ == '__main__':
    main()