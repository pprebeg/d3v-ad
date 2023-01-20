import numpy

from addir.aircraftcomponent import AircraftComponent, AircraftComponentFromMesh
from geometry_extend import GeometryExtension
import numpy as np
import scipy as sp
import openmesh as om
import warnings
import os
def aerop(xg, xd, zg, zd):

    # FUNKCIJA ZA RACUNANJE GEOMETRIJSKIH VELICINA AEROPROFILA

    # Aproksimacija srednje linije polinomom:

    xA = xg
    xA = np.append(xA, xd)
    zA = zg
    zA = np.append(zA, zd)

    p5 = np.polyfit(xA, zA, 5)

    # Aproksimacija profila racionalnom funkcijom:

    def aprf(x, a4, a3, a2, b4, b3, b2, b1, b0):
        a0 = 0
        a1 = -(a4 + a3 + a2 + a0)
        y = (a4 * x ** 4 + a3 * x ** 3 + a2 * x ** 2 + a1 * x + a0) / \
            (b4 * x ** 4 + b3 * x ** 3 + b2 * x ** 2 + b1 * x + b0)
        return y

    x = np.arange(0, 1 + 0.001, 0.001)

    # Aproksimacija gornjake:

    koef = sp.optimize.curve_fit(aprf, xg, zg)
    popt = koef[0]
    koefg = np.array([popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7]])

    yg = aprf(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])

    # Aproksimacija donjake:

    koef = sp.optimize.curve_fit(aprf, xd, zd)
    popt = koef[0]
    koefd = np.array([popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7]])

    yd = aprf(x, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7])

    # Maksimalna debljina:

    t = yg - yd
    dtdx = np.gradient(t, 0.001)

    r = len(dtdx) - 1
    for i in range(r):
        test = dtdx[i] * dtdx[i + 1]
        if test < 0:
            xt = (x[i] + x[i + 1]) / 2
            tc = (t[i] + t[i + 1]) / 2
            break

    # Maksimalna zakrivljenost:

    f = np.polyval(p5, x)

    dfdx = np.gradient(f, 0.001)

    r = len(dfdx) - 1
    for i in range(r):
        test = dfdx[i] * dfdx[i + 1]
        if test < 0:
            xf = (x[i] + x[i + 1]) / 2
            fc = (f[i] + f[i + 1]) / 2
            break

    return p5, koefg, koefd, xt, tc, xf, fc

def bezlen(t, P):

    # FUNKCIJA ZA RACUNANJE DULJINE BEZIEROVE KRIVULJE DRUGOG REDA

    [xi0, xi1, xi2, yi0, yi1, yi2, zi0, zi1, zi2] = P

    if [xi0, yi0, zi0] == [xi1, yi1, zi1] and [xi1, yi1, zi1] == [xi2, yi2, zi2]:
        zeta = 0
        return zeta

    else:
        axi = 2 * (xi0 - 2 * xi1 + xi2)
        ayi = 2 * (yi0 - 2 * yi1 + yi2)
        azi = 2 * (zi0 - 2 * zi1 + zi2)

        bxi = 2 * (-xi0 + xi1)
        byi = 2 * (-yi0 + yi1)
        bzi = 2 * (-zi0 + zi1)

        ai = axi ** 2 + ayi ** 2 + azi ** 2
        bi = 2 * (axi * bxi + ayi * byi + azi * bzi)
        ci = bxi ** 2 + byi ** 2 + bzi ** 2

        warnings.simplefilter('error', RuntimeWarning)
        try:
            zeta = ((ci - (bi ** 2) / (4 * ai)) / np.sqrt(ai)) * \
                   (np.log(np.abs((np.sqrt((t + bi / (2 * ai)) ** 2 + (ci - (bi ** 2) / (4 * ai)) / ai) + t + bi / (2 * ai)) /
                                  (np.sqrt((t + bi / (2 * ai)) ** 2 + (ci - (bi ** 2) / (4 * ai)) / ai) - t - bi / (2 * ai)))) / 4 +
                   (t + bi / (2 * ai)) * np.sqrt((t + bi / (2 * ai)) ** 2 + (ci - (bi ** 2) / (4 * ai)) / ai) / (2 * (ci - (bi ** 2) / (4 * ai)) / ai)
                    - np.log(np.abs((np.sqrt(ci / ai) + bi / (2 * ai)) / (np.sqrt(ci / ai) - bi / (2 * ai)))) / 4
                    - bi * np.sqrt(ci / ai) / (4 * (ci - (bi ** 2) / (4 * ai))))
        except RuntimeWarning:
            zeta = (ai*(t + bi/(2*ai))*np.sqrt((t + bi / (2 * ai)) ** 2 + (ci - (bi ** 2) / (4 * ai)) / ai)/2
                    -bi*np.sqrt(ci/ai)/4)/np.sqrt(ai)

        return zeta
    
def wsurf(s, t, v, xg, xd, zg, zd, ws):

    # FUNKCIJA ZA ODREDIVANJE KOORDINATA KRILA

    import numpy as np

    [HF, cR, H, h, W, w, r1, x22, y22, z22, l, r3, hix51, hiz51, hix52, HW,
     kc12, c22, kc22, c32, kc32, c42, kc42, cT, kc52, a22, ka22, a32, ka32, a42, ka42, a52, ka52] = \
    np.array([v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15],
              v[16], v[17], v[18], v[19], v[20], v[21], v[22], v[23], v[24], v[25], v[26], v[27], v[28], v[29], v[30], v[31], v[32]])

    # Definirane aeroprofila:

    ak = aerop(xg, xd, zg, zd)
    Kc = ak[0]
    Kg = ak[1]
    Kd = ak[2]

    # SEKCIJA 1

    # Krivulja vodilja:

    AB1 = np.array([[1, 0, 0], [0, 0, 1], [0, -1, 1]])
    bBx1 = np.array([HF, H + cR, 0])
    bBy1 = np.array([0, W, r1])
    bBz1 = np.array([0, 0, 0])

    x1kt = np.linalg.solve(AB1, bBx1)
    y1kt = np.linalg.solve(AB1, bBy1)
    z1kt = np.linalg.solve(AB1, bBz1)

    x1 = x1kt[0] * (1 - t) ** 2 + 2 * x1kt[1] * t * (1 - t) + x1kt[2] * t ** 2
    y1 = y1kt[0] * (1 - t) ** 2 + 2 * y1kt[1] * t * (1 - t) + y1kt[2] * t ** 2
    z1 = z1kt[0] * (1 - t) ** 2 + 2 * z1kt[1] * t * (1 - t) + z1kt[2] * t ** 2

    # Distribucija duljine tetive:

    Ac1 = np.array([[0, 0, 0, 1],
                    [w ** 3, w ** 2, w, 1],
                    [W ** 3, W ** 2, W, 1],
                    [3 * W ** 2, 2 * W, 1, 0]])
    bc1 = np.array([H + cR, h + cR, cR, kc12])

    zeta = y1

    Kc1 = np.linalg.solve(Ac1, bc1)
    cG = Kc1[0] * zeta ** 3 + Kc1[1] * zeta ** 2 + Kc1[2] * zeta + Kc1[3]
    ct1 = cG + x1 - (H + cR)

    # Distribucija uvijanja:

    alfa1 = 0

    # Distribucija dihedrala:

    theta1 = 0

    # Krivulja aeroprofila:

    if ws == 'c':
        # Srednja linija
        xa1__ = ct1 * (s - 1)
        za1__ = ct1 * (Kc[0] * s ** 5 + Kc[1] * s ** 4 + Kc[2] * s ** 3 + Kc[3] * s ** 2 + Kc[4] * s + Kc[5])
    elif ws == 'g':
        # Gornjaka
        xa1__ = ct1 * (s - 1)
        za1__ = ct1 * (Kg[0] * s ** 4 + Kg[1] * s ** 3 + Kg[2] * s ** 2 - (Kg[0] + Kg[1] + Kg[2]) * s) / \
                (Kg[3] * s ** 4 + Kg[4] * s ** 3 + Kg[5] * s ** 2 + Kg[6] * s + Kg[7])
    elif ws == 'd':
        # Donjaka
        xa1__ = ct1 * (s - 1)
        za1__ = ct1 * (Kd[0] * s ** 4 + Kd[1] * s ** 3 + Kd[2] * s ** 2 - (Kd[0] + Kd[1] + Kd[2]) * s) / \
                (Kd[3] * s ** 4 + Kd[4] * s ** 3 + Kd[5] * s ** 2 + Kd[6] * s + Kd[7])

    # Transformacija krivulje aeroprofila:

    xa1 = xa1__ * np.cos(alfa1) + za1__ * np.sin(alfa1)
    ya1 = xa1__ * np.sin(alfa1) * np.sin(theta1) - za1__ * np.cos(alfa1) * np.sin(theta1)
    za1 = -xa1__ * np.sin(alfa1) * np.cos(theta1) + za1__ * np.cos(alfa1) * np.cos(theta1)

    # Definiranje plohe:

    Sa1 = np.array([x1 + xa1, y1 + ya1, z1 + za1])

    # SEKCIJA 2

    # Krivulja vodilja:

    AB2 = np.array([[1, 0], [-1, 1]])
    bBx2 = np.array([x1kt[2], -x1kt[1] + x1kt[2]])
    bBy2 = np.array([y1kt[2], -y1kt[1] + y1kt[2]])
    bBz2 = np.array([z1kt[2], -z1kt[1] + z1kt[2]])

    x2kt = np.linalg.solve(AB2, bBx2)
    y2kt = np.linalg.solve(AB2, bBy2)
    z2kt = np.linalg.solve(AB2, bBz2)

    x2 = x2kt[0] * (1 - t) ** 2 + 2 * x2kt[1] * t * (1 - t) + x22 * t ** 2
    y2 = y2kt[0] * (1 - t) ** 2 + 2 * y2kt[1] * t * (1 - t) + y22 * t ** 2
    z2 = z2kt[0] * (1 - t) ** 2 + 2 * z2kt[1] * t * (1 - t) + z22 * t ** 2

    # Distribucija duljine tetive:

    P2 = [x2kt[0], x2kt[1], x22, y2kt[0], y2kt[1], y22, z2kt[0], z2kt[1], z22]
    zeta2 = bezlen(1, P2)

    Ac2 = np.array([[W ** 3, W ** 2, W, 1],
                    [(W + zeta2) ** 3, (W + zeta2) ** 2, W + zeta2, 1],
                    [3 * W ** 2, 2 * W, 1, 0],
                    [3 * (W + zeta2) ** 2, 2 * (W + zeta2), 1, 0]])
    bc2 = np.array([cR, c22, kc12, kc22])

    zeta = W + bezlen(t, P2)

    Kc2 = np.linalg.solve(Ac2, bc2)
    ct2 = Kc2[0] * zeta ** 3 + Kc2[1] * zeta ** 2 + Kc2[2] * zeta + Kc2[3]

    # Distribucija uvijanja:

    Aa2 = np.array([[W ** 3, W ** 2, W, 1],
                    [(W + zeta2) ** 3, (W + zeta2) ** 2, W + zeta2, 1],
                    [3 * W ** 2, 2 * W, 1, 0],
                    [3 * (W + zeta2) ** 2, 2 * (W + zeta2), 1, 0]])
    ba2 = np.array([0, a22, 0, ka22])

    Ka2 = np.linalg.solve(Aa2, ba2)
    alfa2 = Ka2[0] * zeta ** 3 + Ka2[1] * zeta ** 2 + Ka2[2] * zeta + Ka2[3]

    # Distribucija dihedrala:

    ay2 = 2 * (y2kt[0] - 2 * y2kt[1] + y22)
    az2 = 2 * (z2kt[0] - 2 * z2kt[1] + z22)
    by2 = 2 * (-y2kt[0] + y2kt[1])
    bz2 = 2 * (-z2kt[0] + z2kt[1])
    theta2 = np.arctan((az2*t + bz2)/(ay2*t + by2))

    # Krivulja aeroprofila:

    if ws == 'c':
        # Srednja linija
        xa2__ = ct2 * (s - 1)
        za2__ = ct2 * (Kc[0] * s ** 5 + Kc[1] * s ** 4 + Kc[2] * s ** 3 + Kc[3] * s ** 2 + Kc[4] * s + Kc[5])
    elif ws == 'g':
        # Gornjaka
        xa2__ = ct2 * (s - 1)
        za2__ = ct2 * (Kg[0] * s ** 4 + Kg[1] * s ** 3 + Kg[2] * s ** 2 - (Kg[0] + Kg[1] + Kg[2]) * s) / \
                (Kg[3] * s ** 4 + Kg[4] * s ** 3 + Kg[5] * s ** 2 + Kg[6] * s + Kg[7])
    elif ws == 'd':
        # Donjaka
        xa2__ = ct2 * (s - 1)
        za2__ = ct2 * (Kd[0] * s ** 4 + Kd[1] * s ** 3 + Kd[2] * s ** 2 - (Kd[0] + Kd[1] + Kd[2]) * s) / \
                (Kd[3] * s ** 4 + Kd[4] * s ** 3 + Kd[5] * s ** 2 + Kd[6] * s + Kd[7])

    # Transformacija krivulje aeroprofila:

    xa2 = xa2__ * np.cos(alfa2) + za2__ * np.sin(alfa2)
    ya2 = xa2__ * np.sin(alfa2) * np.sin(theta2) - za2__ * np.cos(alfa2) * np.sin(theta2)
    za2 = -xa2__ * np.sin(alfa2) * np.cos(theta2) + za2__ * np.cos(alfa2) * np.cos(theta2)

    # Definiranje plohe:

    Sa2 = np.array([x2 + xa2, y2 + ya2, z2 + za2])

    # SEKCIJA 3

    # Krivulja vodilja:

    xt = (x22 - x2kt[1]) / np.sqrt((x22 - x2kt[1]) ** 2 + (y22 - y2kt[1]) ** 2 + (z22 - z2kt[1]) ** 2)
    yt = (y22 - y2kt[1]) / np.sqrt((x22 - x2kt[1]) ** 2 + (y22 - y2kt[1]) ** 2 + (z22 - z2kt[1]) ** 2)
    zt = (z22 - z2kt[1]) / np.sqrt((x22 - x2kt[1]) ** 2 + (y22 - y2kt[1]) ** 2 + (z22 - z2kt[1]) ** 2)

    AB3 = np.array([[1, 0, 0], [-1, 1, 0], [0, 0, 1]])
    bBx3 = np.array([x22, (l - r3)*xt, x22 + l*xt])
    bBy3 = np.array([y22, (l - r3)*yt, y22 + l*yt])
    bBz3 = np.array([z22, (l - r3)*zt, z22 + l*zt])

    x3kt = np.linalg.solve(AB3, bBx3)
    y3kt = np.linalg.solve(AB3, bBy3)
    z3kt = np.linalg.solve(AB3, bBz3)

    x3 = x3kt[0] * (1 - t) ** 2 + 2 * x3kt[1] * t * (1 - t) + x3kt[2] * t ** 2
    y3 = y3kt[0] * (1 - t) ** 2 + 2 * y3kt[1] * t * (1 - t) + y3kt[2] * t ** 2
    z3 = z3kt[0] * (1 - t) ** 2 + 2 * z3kt[1] * t * (1 - t) + z3kt[2] * t ** 2

    # Distribucija duljine tetive:

    P3 = [x3kt[0], x3kt[1], x3kt[2], y3kt[0], y3kt[1], y3kt[2], z3kt[0], z3kt[1], z3kt[2]]

    Ac3 = np.array([[(W + zeta2) ** 3, (W + zeta2) ** 2, W + zeta2, 1],
                    [(W + zeta2 + l) ** 3, (W + zeta2 + l) ** 2, W + zeta2 + l, 1],
                    [3 * (W + zeta2) ** 2, 2 * (W + zeta2), 1, 0],
                    [3 * (W + zeta2 + l) ** 2, 2 * (W + zeta2 + l), 1, 0]])
    bc3 = np.array([c22, c32, kc22, kc32])

    zeta = W + zeta2 + bezlen(t, P3)

    Kc3 = np.linalg.solve(Ac3, bc3)
    ct3 = Kc3[0] * zeta ** 3 + Kc3[1] * zeta ** 2 + Kc3[2] * zeta + Kc3[3]

    # Distribucija uvijanja:

    Aa3 = np.array([[(W + zeta2) ** 3, (W + zeta2) ** 2, W + zeta2, 1],
                    [(W + zeta2 + l) ** 3, (W + zeta2 + l) ** 2, W + zeta2 + l, 1],
                    [3 * (W + zeta2) ** 2, 2 * (W + zeta2), 1, 0],
                    [3 * (W + zeta2 + l) ** 2, 2 * (W + zeta2 + l), 1, 0]])
    ba3 = np.array([a22, a32, ka22, ka32])

    Ka3 = np.linalg.solve(Aa3, ba3)
    alfa3 = Ka3[0] * zeta ** 3 + Ka3[1] * zeta ** 2 + Ka3[2] * zeta + Ka3[3]

    # Distribucija dihedrala:

    ay3 = 2 * (y3kt[0] - 2 * y3kt[1] + y3kt[2])
    az3 = 2 * (z3kt[0] - 2 * z3kt[1] + z3kt[2])
    by3 = 2 * (-y3kt[0] + y3kt[1])
    bz3 = 2 * (-z3kt[0] + z3kt[1])
    theta3 = np.arctan((az3*t + bz3)/(ay3*t + by3))

    # Krivulja aeroprofila:

    if ws == 'c':
        # Srednja linija
        xa3__ = ct3 * (s - 1)
        za3__ = ct3 * (Kc[0] * s ** 5 + Kc[1] * s ** 4 + Kc[2] * s ** 3 + Kc[3] * s ** 2 + Kc[4] * s + Kc[5])
    elif ws == 'g':
        # Gornjaka
        xa3__ = ct3 * (s - 1)
        za3__ = ct3 * (Kg[0] * s ** 4 + Kg[1] * s ** 3 + Kg[2] * s ** 2 - (Kg[0] + Kg[1] + Kg[2]) * s) / \
                (Kg[3] * s ** 4 + Kg[4] * s ** 3 + Kg[5] * s ** 2 + Kg[6] * s + Kg[7])
    elif ws == 'd':
        # Donjaka
        xa3__ = ct3 * (s - 1)
        za3__ = ct3 * (Kd[0] * s ** 4 + Kd[1] * s ** 3 + Kd[2] * s ** 2 - (Kd[0] + Kd[1] + Kd[2]) * s) / \
                (Kd[3] * s ** 4 + Kd[4] * s ** 3 + Kd[5] * s ** 2 + Kd[6] * s + Kd[7])

    # Transformacija krivulje aeroprofila:

    xa3 = xa3__ * np.cos(alfa3) + za3__ * np.sin(alfa3)
    ya3 = xa3__ * np.sin(alfa3) * np.sin(theta3) - za3__ * np.cos(alfa3) * np.sin(theta3)
    za3 = -xa3__ * np.sin(alfa3) * np.cos(theta3) + za3__ * np.cos(alfa3) * np.cos(theta3)

    # Definiranje plohe:

    Sa3 = np.array([x3 + xa3, y3 + ya3, z3 + za3])

    # SEKCIJA 4

    # Krivulja vodilja:

    AB4 = np.array([[1, 0, 0], [-1, 1, 0], [0, -1, 1]])
    bBx4 = np.array([x3kt[2], -x3kt[1] + x3kt[2], hix51])
    bBy4 = np.array([y3kt[2], -y3kt[1] + y3kt[2], 10**(-10)])
    bBz4 = np.array([z3kt[2], -z3kt[1] + z3kt[2], hiz51])

    x4kt = np.linalg.solve(AB4, bBx4)
    y4kt = np.linalg.solve(AB4, bBy4)
    z4kt = np.linalg.solve(AB4, bBz4)

    x4 = x4kt[0] * (1 - t) ** 2 + 2 * x4kt[1] * t * (1 - t) + x4kt[2] * t ** 2
    y4 = y4kt[0] * (1 - t) ** 2 + 2 * y4kt[1] * t * (1 - t) + y4kt[2] * t ** 2
    z4 = z4kt[0] * (1 - t) ** 2 + 2 * z4kt[1] * t * (1 - t) + z4kt[2] * t ** 2

    # Distribucija duljine tetive:

    P4 = [x4kt[0], x4kt[1], x4kt[2], y4kt[0], y4kt[1], y4kt[2], z4kt[0], z4kt[1], z4kt[2]]
    zeta4 = bezlen(1, P4)

    Ac4 = np.array([[(W + zeta2 + l) ** 3, (W + zeta2 + l) ** 2, W + zeta2 + l, 1],
                    [(W + zeta2 + l + zeta4) ** 3, (W + zeta2 + l + zeta4) ** 2, W + zeta2 + l + zeta4, 1],
                    [3 * (W + zeta2 + l) ** 2, 2 * (W + zeta2 + l), 1, 0],
                    [3 * (W + zeta2 + l + zeta4) ** 2, 2 * (W + zeta2 + l + zeta4), 1, 0]])
    bc4 = np.array([c32, c42, kc32, kc42])

    zeta = W + zeta2 + l + bezlen(t, P4)

    Kc4 = np.linalg.solve(Ac4, bc4)
    ct4 = Kc4[0] * zeta ** 3 + Kc4[1] * zeta ** 2 + Kc4[2] * zeta + Kc4[3]

    # Distribucija uvijanja:

    Aa4 = np.array([[(W + zeta2 + l) ** 3, (W + zeta2 + l) ** 2, W + zeta2 + l, 1],
                    [(W + zeta2 + l + zeta4) ** 3, (W + zeta2 + l + zeta4) ** 2, W + zeta2 + l + zeta4, 1],
                    [3 * (W + zeta2 + l) ** 2, 2 * (W + zeta2 + l), 1, 0],
                    [3 * (W + zeta2 + l + zeta4) ** 2, 2 * (W + zeta2 + l + zeta4), 1, 0]])
    ba4 = np.array([a32, a42, ka32, ka42])

    Ka4 = np.linalg.solve(Aa4, ba4)
    alfa4 = Ka4[0] * zeta ** 3 + Ka4[1] * zeta ** 2 + Ka4[2] * zeta + Ka4[3]

    # Distribucija dihedrala:

    ay4 = 2 * (y4kt[0] - 2 * y4kt[1] + y4kt[2])
    az4 = 2 * (z4kt[0] - 2 * z4kt[1] + z4kt[2])
    by4 = 2 * (-y4kt[0] + y4kt[1])
    bz4 = 2 * (-z4kt[0] + z4kt[1])
    theta4 = np.arctan((az4*t + bz4)/(ay4*t + by4))

    # Krivulja aeroprofila:

    if ws == 'c':
        # Srednja linija
        xa4__ = ct4 * (s - 1)
        za4__ = ct4 * (Kc[0] * s ** 5 + Kc[1] * s ** 4 + Kc[2] * s ** 3 + Kc[3] * s ** 2 + Kc[4] * s + Kc[5])
    elif ws == 'g':
        # Gornjaka
        xa4__ = ct4 * (s - 1)
        za4__ = ct4 * (Kg[0] * s ** 4 + Kg[1] * s ** 3 + Kg[2] * s ** 2 - (Kg[0] + Kg[1] + Kg[2]) * s) / \
                (Kg[3] * s ** 4 + Kg[4] * s ** 3 + Kg[5] * s ** 2 + Kg[6] * s + Kg[7])
    elif ws == 'd':
        # Donjaka
        xa4__ = ct4 * (s - 1)
        za4__ = ct4 * (Kd[0] * s ** 4 + Kd[1] * s ** 3 + Kd[2] * s ** 2 - (Kd[0] + Kd[1] + Kd[2]) * s) / \
                (Kd[3] * s ** 4 + Kd[4] * s ** 3 + Kd[5] * s ** 2 + Kd[6] * s + Kd[7])

    # Transformacija krivulje aeroprofila:

    xa4 = xa4__ * np.cos(alfa4) + za4__ * np.sin(alfa4)
    ya4 = xa4__ * np.sin(alfa4) * np.sin(theta4) - za4__ * np.cos(alfa4) * np.sin(theta4)
    za4 = -xa4__ * np.sin(alfa4) * np.cos(theta4) + za4__ * np.cos(alfa4) * np.cos(theta4)

    # Definiranje plohe:

    Sa4 = np.array([x4 + xa4, y4 + ya4, z4 + za4])

    # SEKCIJA 5

    # Krivulja vodilja:

    AB5 = np.array([[1, 0, 0], [-1, 1, 0], [0, 0, 1]])
    bBx5 = np.array([x4kt[2], -x4kt[1] + x4kt[2], x4kt[2] + hix52])
    bBy5 = np.array([y4kt[2], -y4kt[1] + y4kt[2], y4kt[2]])
    bBz5 = np.array([z4kt[2], -z4kt[1] + z4kt[2], z4kt[2] + HW])

    x5kt = np.linalg.solve(AB5, bBx5)
    y5kt = np.linalg.solve(AB5, bBy5)
    z5kt = np.linalg.solve(AB5, bBz5)

    x5 = x5kt[0] * (1 - t) ** 2 + 2 * x5kt[1] * t * (1 - t) + x5kt[2] * t ** 2
    y5 = y5kt[0] * (1 - t) ** 2 + 2 * y5kt[1] * t * (1 - t) + y5kt[2] * t ** 2
    z5 = z5kt[0] * (1 - t) ** 2 + 2 * z5kt[1] * t * (1 - t) + z5kt[2] * t ** 2

    # Distribucija duljine tetive:

    P5 = [x5kt[0], x5kt[1], x5kt[2], y5kt[0], y5kt[1], y5kt[2], z5kt[0], z5kt[1], z5kt[2]]
    zeta5 = bezlen(1, P5)

    Ac5 = np.array([[(W + zeta2 + l + zeta4) ** 3, (W + zeta2 + l + zeta4) ** 2, W + zeta2 + l + zeta4, 1],
                    [(W + zeta2 + l + zeta4 + zeta5) ** 3, (W + zeta2 + l + zeta4 + zeta5) ** 2, W + zeta2 + l + zeta4 + zeta5, 1],
                    [3 * (W + zeta2 + l + zeta4) ** 2, 2 * (W + zeta2 + l + zeta4), 1, 0],
                    [3 * (W + zeta2 + l + zeta4 + zeta5) ** 2, 2 * (W + zeta2 + l + zeta4 + zeta5), 1, 0]])
    bc5 = np.array([c42, cT, kc42, kc52])

    zeta = W + zeta2 + l + zeta4 + bezlen(t, P5)

    Kc5 = np.linalg.solve(Ac5, bc5)
    ct5 = Kc5[0] * zeta ** 3 + Kc5[1] * zeta ** 2 + Kc5[2] * zeta + Kc5[3]

    # Distribucija uvijanja:

    Aa5 = np.array([[(W + zeta2 + l + zeta4) ** 3, (W + zeta2 + l + zeta4) ** 2, W + zeta2 + l + zeta4, 1],
                    [(W + zeta2 + l + zeta4 + zeta5) ** 3, (W + zeta2 + l + zeta4 + zeta5) ** 2,
                     W + zeta2 + l + zeta4 + zeta5, 1],
                    [3 * (W + zeta2 + l + zeta4) ** 2, 2 * (W + zeta2 + l + zeta4), 1, 0],
                    [3 * (W + zeta2 + l + zeta4 + zeta5) ** 2, 2 * (W + zeta2 + l + zeta4 + zeta5), 1, 0]])
    ba5 = np.array([a42, a52, ka42, ka52])

    Ka5 = np.linalg.solve(Aa5, ba5)
    alfa5 = Ka5[0] * zeta ** 3 + Ka5[1] * zeta ** 2 + Ka5[2] * zeta + Ka5[3]

    # Distribucija dihedrala:

    theta5 = np.pi/2

    # Krivulja aeroprofila:

    if ws == 'c':
        # Srednja linija
        xa5__ = ct5 * (s - 1)
        za5__ = ct5 * (Kc[0] * s ** 5 + Kc[1] * s ** 4 + Kc[2] * s ** 3 + Kc[3] * s ** 2 + Kc[4] * s + Kc[5])
    elif ws == 'g':
        # Gornjaka
        xa5__ = ct5 * (s - 1)
        za5__ = ct5 * (Kg[0] * s ** 4 + Kg[1] * s ** 3 + Kg[2] * s ** 2 - (Kg[0] + Kg[1] + Kg[2]) * s) / \
                (Kg[3] * s ** 4 + Kg[4] * s ** 3 + Kg[5] * s ** 2 + Kg[6] * s + Kg[7])
    elif ws == 'd':
        # Donjaka
        xa5__ = ct5 * (s - 1)
        za5__ = ct5 * (Kd[0] * s ** 4 + Kd[1] * s ** 3 + Kd[2] * s ** 2 - (Kd[0] + Kd[1] + Kd[2]) * s) / \
                (Kd[3] * s ** 4 + Kd[4] * s ** 3 + Kd[5] * s ** 2 + Kd[6] * s + Kd[7])

    # Transformacija krivulje aeroprofila:

    xa5 = xa5__ * np.cos(alfa5) + za5__ * np.sin(alfa5)
    ya5 = xa5__ * np.sin(alfa5) * np.sin(theta5) - za5__ * np.cos(alfa5) * np.sin(theta5)
    za5 = -xa5__ * np.sin(alfa5) * np.cos(theta5) + za5__ * np.cos(alfa5) * np.cos(theta5)

    # Definiranje plohe:

    Sa5 = np.array([x5 + xa5, y5 + ya5, z5 + za5])

    return Sa1, Sa2, Sa3, Sa4, Sa5

class BraMoFlyingWing(AircraftComponent):
    def __init__(self, fileName, name=''):
        super().__init__(name)
        self._filename = fileName
        self._v,self._afoil_xy = self.read_file()
        #self.regenerate_component()
        self.regenerate_subcomponents()

    @property
    def filename(self):
        return self._filename

    def read_file(self):
        v=[0.0]*33
        v[0] = 5  # HF
        v[1] = 3.6  # cR
        v[2] = 0.8  # H
        v[3] = 0.4  # h
        v[4] = 2  # W
        v[5] = 1  # w
        v[6] = 0.6  # r1
        v[7] = 5  # x22
        v[8] = 4.46  # y22
        v[9] = 0.2  # z22
        v[10] = 6  # l
        v[11] = 0.6  # r3
        v[12] = 0.24  # hix51
        v[13] = 0.2  # hiz51
        v[14] = 0.3  # hix52
        v[15] = 0.5  # HW
        v[16] = -0.4  # kc12
        v[17] = 2.6  # c22
        v[18] = -0.5  # kc22
        v[19] = 0.9  # c32 0.9
        v[20] = -0.5  # kc32
        v[21] = 0.5  # c42
        v[22] = -0.6  # kc42
        v[23] = 0.02  # cT
        v[24] = -1  # kc52
        v[25] = (0) * np.pi / 180  # a22
        v[26] = -0.1 * 0  # ka22
        v[27] = (0) * np.pi / 180  # a32
        v[28] = -0.1 * 0  # ka32
        v[29] = (0) * np.pi / 180  # a42
        v[30] = -0.1 * 0  # ka42
        v[31] = (0) * np.pi / 180  # a52
        v[32] = -0.1 * 0  # ka52
        afoil_xy = np.loadtxt(self._filename)

        return v,afoil_xy


    def get_mesh_data(self):
        v = self._v
        MH60 = self._afoil_xy
        xg = MH60[0:57, 0]
        xd = MH60[57:, 0]
        zg = MH60[0:57, 1]
        zd = MH60[57:, 1]

        s = np.arange(0, 1 + 0.01, 0.01)
        t = np.arange(0, 1 + 0.01, 0.01)
        [S, T] = np.meshgrid(s, t)
        Sc = wsurf(S, T, v, xg, xd, zg, zd, 'c')
        Sg = wsurf(S, T, v, xg, xd, zg, zd, 'g')
        Sd = wsurf(S, T, v, xg, xd, zg, zd, 'd')
        return Sc,Sg,Sd

    def regenerate_component(self):
        Sc,Sg,Sd = self.get_mesh_data()
        shape = np.shape(Sc)
        nsec= shape[0]
        nc = shape[1]
        nl = shape[2]
        points = []
        fvi = []
        for sec in range(nsec):
            points.extend(np.transpose(np.reshape(Sc[sec], (nc, nl * nl))))
            for iwl in range(nl-1):
                for ip in range(nl-1-1):
                    ip_wl_0_p_0 = sec*(nl*nl)+iwl+ip*nl
                    ip_wl_0_p_1 = ip_wl_0_p_0+1
                    ip_wl_1_p_0 = sec*(nl*nl)+(ip+1)*nl+iwl
                    ip_wl_1_p_1 = ip_wl_1_p_0+1
                    fvi.append([ip_wl_0_p_0,ip_wl_0_p_1,ip_wl_1_p_0])
                    fvi.append([ip_wl_1_p_0, ip_wl_0_p_1, ip_wl_1_p_1])
        #points = np.reshape(points,(nsec*nl*nl,nc))
        self.mesh = om.TriMesh(points, fvi)

    def regenerate_subcomponents(self):
        S = self.get_mesh_data()
        names = ['Mid plane', 'Upper surface','Lower surface']
        for ipart in range(len(S)):
            shape = np.shape(S[ipart])
            nsec= shape[0]
            nc = shape[1]
            nl = shape[2]
            for sec in range(nsec):
                fvi = []
                points = np.transpose(np.reshape(S[ipart][sec], (nc, nl * nl)))
                for iwl in range(nl-1):
                    for ip in range(nl-1):
                        ip_wl_0_p_0 = iwl+ip*nl
                        ip_wl_0_p_1 = ip_wl_0_p_0+1
                        ip_wl_1_p_0 = (ip+1)*nl+iwl
                        ip_wl_1_p_1 = ip_wl_1_p_0+1
                        fvi.append([ip_wl_0_p_0,ip_wl_0_p_1,ip_wl_1_p_0])
                        fvi.append([ip_wl_1_p_0, ip_wl_0_p_1, ip_wl_1_p_1])
                sg = GeometryExtension(names[ipart]+ ' Section '+str(sec+1))
                sg.mesh=om.TriMesh(points, fvi)
                sg.emit_geometry_built()

    def exportGeometry(self, fileName):
        AircraftComponentFromMesh.export_component(fileName, self)