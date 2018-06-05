import numpy
import json
import glob

def main():
    folder = "random2"
    data = {
    "steps": 2000,
    "seed": 1886615,
    "dr": 1.0,
    "gp": 64,
    "dt": 0.05,
    "L": 63.0,
    "N": 1000000,
    "vt": 1.0,
    "vd": 5.0,
    "field": {
        "Bx": 0.0,
        "By": 0.0,
        "Bz": 0.1
        },
    "Efield": {
        "Ex": 0.0,
        "Ey": 0.05,
        "Ez": 0.0
        }
    }
    outfile = "../{}.json".format(folder)
    data["Directory"] = "/home/daniel/Desktop/random_particles/{}".format(folder)
    with open(outfile, "w") as file:
        json.dump(data, file)

if __name__ == '__main__':
    main()