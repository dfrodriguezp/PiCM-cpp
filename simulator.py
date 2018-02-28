import time
import numpy
import particles
import os
import functions
import click
import json
from tqdm import tqdm
from matplotlib import pyplot
import getpass

@click.command()
@click.option("--system", type=str, default="TSI")
@click.option("--test", type=str, default="1")
@click.option("--c", type=int, prompt="C", default=2)
@click.option("--np", type=int, prompt="np", default=5, help="Cubic root of number of particles")
@click.option("--dt", type=float, prompt="dt", default=0.1)
# @click.option("--dr", type=float, prompt="dr", default=1.0)
@click.option("--l", type=float, prompt="L", default=1.0)
@click.option("--steps", type=int, prompt="steps", default=10)
@click.option("--seed", type=int, prompt="seed", default=6969)
def main(system, test, c, np, dt, l, steps, seed):
    GP = c + 1
    NP = np ** 2
    L = l
    dr = L / c

    vt = L * 0.0015626
    # vd = L * 0.03125
    vd = 0.05
    # B = numpy.array([0.0, 1.0])

    if system == "TSI":
        parts = particles.two_stream(NP, L, dr, vt, vd, seed)

    # directory = "/home/%s/Desktop/%s/Test_%s" %(getpass.getuser(), system, test)
    # parameters = "/home/%s/Desktop/%s/parameters" %(getpass.getuser(), system)
    # folders = ["/phase_space", "/velocities", "/rho", "/phi", "/Efield", "/energies"]

    # for f in folders:
    #     os.system("mkdir -p %s" %(directory + f))

    # sample = {
    # "directory":directory,
    # "system":system,
    # "grid points": (GP, GP),
    # "number of particles":NP,
    # "box size":L,
    # "steps":steps
    # }

    # os.system("mkdir -p %s" %parameters)

    # with open("%s/%s_Test_%s.json" %(parameters, system, test), "w") as fp:
    #     json.dump(sample, fp)

    # energies = open("%s/energies/energies.dat" %directory, "w")

    for step in tqdm(range(steps)):
        RHO = functions.density2(parts, GP, dr)
        PHIn, PHI_k, RHO_k = functions.potential(RHO, GP, dr)
        EFIELDn = functions.Efield_GP(PHIn, GP, dr)
        EFIELDp = functions.Efield_P(EFIELDn, parts, dr)
    
        if step == 0:
            parts = functions.rewind(-1, EFIELDp, parts, dt)
    
        parts = functions.leapfrog(EFIELDp, parts, L, dt)
        final_parts = functions.rewind(1, EFIELDp, parts, dt)

        pos = numpy.array([p.pos for p in final_parts])
        print(pos)
    #
    #     phase_space = open("%s/phase_space/step%s.dat" %(directory, step), "w")
    #     velocities = open("%s/velocities/step%s.dat" %(directory, step), "w")
    #     E = open("%s/Efield/step%s.dat" %(directory, step), "w")
    #     phi = open("%s/phi/step%s.dat" %(directory, step), "w")
    #     rho = open("%s/rho/step%s.dat" %(directory, step), "w")
    #
    #     for p in final_parts:
    #         if p.move:
    #             phase_space.write("%s %s" %(p.pos[0], p.vel[0]))
    #             phase_space.write("\n")
    #             velocities.write("%s %s" %(p.vel[0], p.vel[1]))
    #             velocities.write("\n")
    #     E_F = 0
    #     for i in range(GP):
    #         for j in range(GP):
    #             E.write("%s %s %s %s" %(i*dr, j*dr, EFIELDn[0][i][j], EFIELDn[1][i][j]))
    #             E.write("\n")
    #             phi.write("%s %s %s" %(i*dr, j*dr, PHIn[i][j]))
    #             phi.write("\n")
    #             rho.write("%s %s %s" %(i*dr, j*dr, RHO[i][j]))
    #             rho.write("\n")
    #             E_F += RHO[i][j] * PHIn[i][j]
    #     E_F *= 0.5
    #
    #     energies.write("%s %s" %(E_F, step))
    #     energies.write("\n")
    #     phase_space.close()
    #     E.close()
    #     phi.close()
    #     rho.close()
    #     velocities.close()
    # energies.close()

if __name__ == '__main__':
    main()
