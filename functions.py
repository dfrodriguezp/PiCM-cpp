import numpy

def density1(parts, gp, dr, rho_c):
    rho = numpy.zeros(shape=(gp, gp))

    pos = numpy.array([p.pos for p in parts])
    node = numpy.array(pos / dr, dtype=int)

    h = pos - node*dr

    A = rho_c * (dr - h[:, 0]) * (dr - h[:, 1])
    B = rho_c * (dr - h[:, 0]) * h[:, 1]
    C = rho_c * h[:, 0] * (dr - h[:, 1])
    D = rho_c * h[:, 0] * h[:, 1]

    nodes_list = [(i, j) for i in range(gp-1) for j in range(gp-1)]

    for n in nodes_list:
        i, j = n
        where = numpy.all(node == n, axis=1)
        rho[i, j] += numpy.sum(A[where])
        rho[i, j + 1] += numpy.sum(B[where])
        rho[i + 1, j] += numpy.sum(C[where])
        rho[i + 1, j + 1] += numpy.sum(D[where])

    rho[:, -1] = (rho[:, -1] + rho[:, 0]) * 0.5
    rho[-1, :] = (rho[-1, :] + rho[0, :]) * 0.5
    rho[:, 0] = rho[:, -1]
    rho[0, :] = rho[-1, :]

    rho /= (dr * dr)

    return rho

def density2(parts, gp, dr, rho_c):
    rho = numpy.zeros(shape=(gp, gp))

    for ind, p in enumerate(parts):
        i = int(p.pos[0] / dr)
        j = int(p.pos[1] / dr)

        hx = p.pos[0] - (i * dr)
        hy = p.pos[1] - (j * dr)

        # rho_c = p.q / (dr * dr)

        rho[i, j] += rho_c[ind] * (dr - hx) * (dr - hy)
        rho[i, j+1] += rho_c[ind] * (dr - hx) * (hy)
        rho[i+1, j] += rho_c[ind] * (hx) * (dr - hy)
        rho[i+1, j+1] += rho_c[ind] * (hx) * (hy)


    rho[:, -1] = (rho[:, -1] + rho[:, 0]) * 0.5
    rho[-1, :] = (rho[-1, :] + rho[0, :]) * 0.5
    rho[:, 0] = rho[:, -1]
    rho[0, :] = rho[-1, :]

    rho /= (dr * dr)
    
    return rho

def potential(rho, gp, dr):
    rho_k = numpy.fft.fft2(rho)
    phi_k = numpy.copy(rho_k)
    i = 0.0 + 1.0j
    W = numpy.exp(2 * numpy.pi * i / gp)
    Wm = 1.0 + 0.0j
    Wn = 1.0 + 0.0j

    for m in range(gp):
        for n in range(gp):
            denom = 4.0 + 0.0j
            denom -= Wm + 1.0/Wm + Wn + 1.0/Wn
            if denom != 0.0:
                phi_k[m][n] *= dr*dr / denom
            Wn *= W
        Wm *= W

    phi = numpy.fft.ifft2(phi_k)
    phi = numpy.real(phi)

    return phi, phi_k, rho_k

def Efield_GP(phi, gp, dr):
    E = numpy.zeros(shape=(gp, gp, 2))

    for j in range(gp):
        for i in range(gp):
            nxt_i = (i + 1) if (i < gp-1) else 0
            prv_i = (i - 1) if (i > 0) else (gp-1)

            E[i][j][0] = (phi[prv_i][j] - phi[nxt_i][j]) / (2 * dr)

    for i in range(gp):
        for j in range(gp):
            nxt_j = (j + 1) if (j < gp-1) else 0
            prv_j = (j - 1) if (j > 0) else (gp-1)

            E[i][j][1] = (phi[i][prv_j] - phi[i][nxt_j]) / (2 * dr)

    return E

def Efield_P(field, parts, dr):
    E = numpy.zeros(shape=(len(parts), 2))

    index = 0
    for p in parts:
        if p.move:
            i = int(p.pos[0] / dr)
            j = int(p.pos[1] / dr)
            hx = p.pos[0] - (i * dr)
            hy = p.pos[1] - (j * dr)

            E[index][0] = field[i][j][0] * (dr - hx) * (dr - hy) + field[i][j+1][0] * (dr - hx) * (hy) + field[i+1][j][0] * (hx) * (dr - hy) + field[i+1][j+1][0] * (hx) * (hy)
            E[index][1] = field[i][j][1] * (dr - hx) * (dr - hy) + field[i][j+1][1] * (dr - hx) * (hy) + field[i+1][j][1] * (hx) * (dr - hy) + field[i+1][j+1][1] * (hx) * (hy)

        index += 1

    E /= (dr * dr)

    return E

def Boris(E, B, parts, L, dt):
    index = 0
    for p in parts:
        if p.move:
            t = 0.5 * (p.qm) * B * dt
            t_2 = numpy.linalg.norm(t) * numpy.linalg.norm(t)
            s = (2 * t) / (1 + t_2)
            v_minus = p.vel + 0.5 * (p.qm) * E[index] * dt
            v_prime = v_minus + numpy.cross(v_minus, t)
            v_plus = v_minus + numpy.cross(v_prime, s)
            p.vel = v_plus + 0.5 * (p.qm) * E[index] * dt

            p.pos += p.vel * dt

            p.pos = p.pos % L

        index += 1

    return parts

def leapfrog(E, parts, L, dt):
    index = 0
    for p in parts:
        if p.move:
            p.vel += (p.qm) * E[index] * dt
            p.pos += p.vel * dt

            p.pos = p.pos % L

        index += 1
    return parts

# def rewind(direction, E, parts, dt):
#     index = 0
#     for p in parts:
#         if p.move:
#             p.vel += direction * 0.5 * (p.qm) * E[index] * dt
#         index += 1
#     return parts

def rewind(direction, E, B, parts, dt):
    dt = dt * 0.5
    index = 0
    for p in parts:
        if p.move:
            t = 0.5 * (p.qm) * B * dt
            t_2 = numpy.linalg.norm(t) * numpy.linalg.norm(t)
            s = (2 * t) / (1 + t_2)
            v_minus = p.vel + 0.5 * (p.qm) * E[index] * dt
            v_prime = v_minus + numpy.cross(v_minus, t)
            v_plus = v_minus + numpy.cross(v_prime, s)
            p.vel = direction * (v_plus + 0.5 * (p.qm) * E[index] * dt)

        index += 1

    return parts
