import time
import numpy
import particles
import os
import functions
import click
import json
from tqdm import tqdm
import getpass

@click.command()
@click.option("--system", type=str, default="TSI")
# @click.option("--test", type=str, default="1")
@click.option("--bx", type=float, default=0.0)
@click.option("--by", type=float, default=0.0)
@click.option("--gp", type=int, default=4)
@click.option("--np", type=int, default=3)
@click.option("--dt", type=float, default=0.1)
# @click.option("--dr", type=float, prompt="dr", default=1.0)
@click.option("--l", type=float, default=1.0)
@click.option("--p_vt", type=float, default=0.001526)
@click.option("--p_vd", type=float, default=0.03125)
@click.option("--steps", type=int, default=1)
@click.option("--seed", type=int, default=69)
def main(system, bx, by, gp, np, dt, l, p_vt, p_vd, steps, seed):
    GP = gp
    NP = np ** 2
    L = l
    dr = L / (gp - 1)

    vt = L * p_vt
    vd = L * p_vd
    B = numpy.array([bx, by])

    if system == "TSI":
        parts = particles.two_stream(NP, L, dr, vt, vd, seed)

    cosa = "GP_{}_NP_{}_dr_{}_Bx_{}_By_{}_vt_{}_vd_{}".format(GP, NP, dr, bx, by, vt, vd)
    directory = "/home/{}/Desktop/{}/{}".format(getpass.getuser(), system, cosa)
    parameters = "/home/{}/Desktop/{}/parameters".format(getpass.getuser(), system)
    folders = ["/space", "/phase_space_x", "/phase_space_y", "/velocities", "/rho", "/phi", "/Efield", "/energies"]

    for f in folders:
        os.system("mkdir -p {}".format(directory + f))

    sample = {
    "directory": directory,
    "system": system,
    "grid points": (GP, GP),
    "number of particles": NP,
    "box size": L,
    "vt": vt,
    "vd": vd,
    "B": (bx, by),
    "steps": steps
    }

    os.system("mkdir -p {}".format(parameters))

    with open("{}/{}_{}.json".format(parameters, system, cosa), "w") as fp:
        json.dump(sample, fp)

    energies = open("{}/energies/energies.dat".format(directory), "w")
    rho_c = numpy.array([p.q / (dr * dr) for p in parts])
    for step in tqdm(range(steps)):
        RHO = functions.density2(parts, GP, dr, rho_c)
        PHIn, PHI_k, RHO_k = functions.potential(RHO, GP, dr)
        EFIELDn = functions.Efield_GP(PHIn, GP, dr)
        EFIELDp = functions.Efield_P(EFIELDn, parts, dr)
    
        if step == 0:
            parts = functions.rewind(-1, EFIELDp, B, parts, dt)
    
        parts = functions.Boris(EFIELDp, B, parts, L, dt)
        final_parts = functions.rewind(1, EFIELDp, B, parts, dt)
    
        phase_space_x = open("{}/phase_space_x/step{}.dat".format(directory, step), "w")
        phase_space_y = open("{}/phase_space_y/step{}.dat".format(directory, step), "w")
        space = open("{}/space/step{}.dat".format(directory, step), "w")
        velocities = open("{}/velocities/step{}.dat".format(directory, step), "w")
        E = open("{}/Efield/step{}.dat".format(directory, step), "w")
        phi = open("{}/phi/step{}.dat".format(directory, step), "w")
        rho = open("{}/rho/step{}.dat".format(directory, step), "w")
    
        for p in final_parts:
            if p.move:
                phase_space_x.write("{} {}".format(p.pos[0], p.vel[0]))
                phase_space_x.write("\n")
                phase_space_y.write("{} {}".format(p.pos[1], p.vel[1]))
                phase_space_y.write("\n")
                space.write("{} {}".format(p.pos[0], p.pos[1]))
                space.write("\n")
                velocities.write("{} {}".format(p.vel[0], p.vel[1]))
                velocities.write("\n")
        vels = numpy.array([p.vel for p in final_parts])
        E_F = 0
        E_k = 0.5 * NP * numpy.sum(numpy.linalg.norm(vels, axis=1) * numpy.linalg.norm(vels, axis=1))
        for i in range(GP):
            for j in range(GP):
                E.write("{} {} {} {}".format(i*dr, j*dr, EFIELDn[i][j][0], EFIELDn[i][j][1]))
                E.write("\n")
                phi.write("{} {} {}".format(i*dr, j*dr, PHIn[i][j]))
                phi.write("\n")
                rho.write("{} {} {}".format(i*dr, j*dr, RHO[i][j]))
                rho.write("\n")
                E_F += RHO[i][j] * PHIn[i][j]
        E_F *= 0.5
    
        energies.write("{} {} {}".format(E_k, E_F, step))
        energies.write("\n")
        phase_space_x.close()
        phase_space_y.close()
        space.close()
        E.close()
        phi.close()
        rho.close()
        velocities.close()
    energies.close()

if __name__ == '__main__':
    main()
