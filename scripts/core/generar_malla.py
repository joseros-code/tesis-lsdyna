
# generar_malla.py
# ------------------------------------------------------------
# Malla estructurada nucleo-citoplasma.
#
# La celula se construye desde un cubo central cartesiano y luego
# se extiende radialmente hasta la superficie externa. La densidad d
# controla directamente las separaciones del cubo central.
# ------------------------------------------------------------

import argparse
import math
from pathlib import Path

from definir_materiales import (
    args_mat,
    leer_mat,
    mat_block,
)


base = Path(__file__).resolve().parents[2]


def leer_args():
    p = argparse.ArgumentParser(description="Genera la malla nucleo-citoplasma.")
    p.add_argument("--Rc", "--r", dest="Rc", type=float, default=0.005)
    p.add_argument("--eta", type=float, default=0.67)
    p.add_argument("--densidad", "--d", dest="densidad", type=int, default=6)
    p.add_argument("--cx", type=float, default=0.0)
    p.add_argument("--cy", type=float, default=0.0)
    p.add_argument("--cz", type=float, default=0.01)
    p.add_argument(
        "--q-cubo",
        type=float,
        default=0.35,
        help="Tamano relativo del cubo central respecto de Rc.",
    )
    p.add_argument("--salida", "--out", dest="salida", default="")
    args_mat(p)
    return p.parse_args()


args = leer_args()
mat = leer_mat(args)


# Datos principales
pid_cito = 2
pid_nuc = 5
mid_cito = 2
mid_nuc = 5
secid = 2

cx = args.cx
cy = args.cy
cz = args.cz

Rc = args.Rc
eta = args.eta
Rn = Rc * eta
d = args.densidad
q_cubo = min(args.q_cubo, eta / math.sqrt(3.0) * 0.95, 0.55)
a = Rc * q_cubo


def lin(a, b, n):
    vals = []

    if n <= 0:
        return [a]

    for i in range(n + 1):
        t = i / n
        vals.append(a + (b - a) * t)

    return vals


def divisiones_radiales():
    # d controla directamente las separaciones del cubo central.
    # El cascaron usa 2d capas radiales, lo que reproduce el conteo
    # de elementos de la malla original de Carlos: 14*d^3.
    if d < 1:
        raise ValueError("La densidad debe ser mayor o igual a 1.")

    if not (0.0 < eta < 1.0):
        raise ValueError("eta debe estar entre 0 y 1.")

    if not (0.0 < q_cubo < eta):
        raise ValueError("q_cubo debe ser mayor que 0 y menor que eta.")

    capas_shell = 2 * d
    capas_nuc = max(1, int(round(capas_shell * (eta - q_cubo) / (1.0 - q_cubo))))

    if capas_nuc >= capas_shell:
        capas_nuc = capas_shell - 1

    capas_cito = capas_shell - capas_nuc

    return capas_nuc, capas_cito


capas_nuc, capas_cito = divisiones_radiales()


capas_radiales = capas_nuc + capas_cito


def punto_radial(px, py, pz, capa):
    # px, py, pz son coordenadas locales de un punto sobre la frontera
    # del cubo. Cada linea radial conserva esa direccion hasta la esfera.
    r0 = math.sqrt(px * px + py * py + pz * pz)
    dx = px / r0
    dy = py / r0
    dz = pz / r0

    if capa == 0:
        r = r0
    elif capa <= capas_nuc:
        t = capa / capas_nuc
        r = r0 + (Rn - r0) * t
    else:
        t = (capa - capas_nuc) / capas_cito
        r = Rn + (Rc - Rn) * t

    return cx + r * dx, cy + r * dy, cz + r * dz


def nodo_txt(nid, xyz):
    x, y, z = xyz
    return f"{nid:8d}{x:16.8E}{y:16.8E}{z:16.8E}\n"


def elem_txt(eid, pid, conn):
    txt = f"{eid:8d}{pid:8d}"

    for nid in conn:
        txt = txt + f"{nid:8d}"

    return txt + "\n"


def sec_solid():
    return [
        "$\n",
        "$ Seccion solida\n",
        "$\n",
        "*SECTION_SOLID\n",
        "$#   secid    elform       aet\n",
        f"{secid:10d}{1:10d}{0:10d}\n",
    ]


def part(pid, mid, nombre):
    return [
        "$\n",
        f"$ Parte {nombre}\n",
        "$\n",
        "*PART\n",
        f"{nombre}\n",
        "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n",
        f"{pid:10d}{secid:10d}{mid:10d}"
        f"{0:10d}{0:10d}{0:10d}{0:10d}{0:10d}\n",
    ]


if args.salida:
    salida = Path(args.salida)

    if not salida.is_absolute():
        salida = base / salida
else:
    rtag = f"{Rc:.4f}".replace(".", "p")
    etag = f"{eta:.3f}".replace(".", "p")
    salida = base / "generados" / f"Celula_basic_D{d:02d}_Rc{rtag}_eta{etag}.k"


# Nodos y elementos
nodos = {}
ids_por_coord = {}
nid = 1

elems = []
eid = 1


def clave(xyz):
    return tuple(round(v, 12) for v in xyz)


def add_nodo(xyz):
    global nid
    key = clave(xyz)

    if key in ids_por_coord:
        return ids_por_coord[key]

    ids_por_coord[key] = nid
    nodos[nid] = xyz
    nid = nid + 1

    return ids_por_coord[key]


def resta(a, b):
    return a[0] - b[0], a[1] - b[1], a[2] - b[2]


def cruz(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def punto(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def orientar(conn):
    p1 = nodos[conn[0]]
    p2 = nodos[conn[1]]
    p4 = nodos[conn[3]]
    p5 = nodos[conn[4]]
    det = punto(resta(p2, p1), cruz(resta(p4, p1), resta(p5, p1)))

    if det < 0.0:
        conn = [conn[0], conn[3], conn[2], conn[1], conn[4], conn[7], conn[6], conn[5]]

    return conn


def add_elem(pid, conn):
    global eid
    elems.append((eid, pid, orientar(conn)))
    eid = eid + 1


xs = lin(-a, a, 2 * d)
ys = lin(0.0, a, d)
zs = lin(-a, 0.0, d)

# Cubo central.
id_cubo = {}

for i, x in enumerate(xs):
    for j, y in enumerate(ys):
        for k, z in enumerate(zs):
            id_cubo[(i, j, k)] = add_nodo((cx + x, cy + y, cz + z))

for i in range(2 * d):
    for j in range(d):
        for k in range(d):
            conn = [
                id_cubo[(i, j, k)],
                id_cubo[(i + 1, j, k)],
                id_cubo[(i + 1, j + 1, k)],
                id_cubo[(i, j + 1, k)],
                id_cubo[(i, j, k + 1)],
                id_cubo[(i + 1, j, k + 1)],
                id_cubo[(i + 1, j + 1, k + 1)],
                id_cubo[(i, j + 1, k + 1)],
            ]
            add_elem(pid_nuc, conn)


def pid_capa(capa):
    return pid_nuc if capa < capas_nuc else pid_cito


def bloque_x(signo):
    ids = {}

    for l in range(capas_radiales + 1):
        for j, y in enumerate(ys):
            for k, z in enumerate(zs):
                ids[(l, j, k)] = add_nodo(punto_radial(signo * a, y, z, l))

    for l in range(capas_radiales):
        for j in range(d):
            for k in range(d):
                conn = [
                    ids[(l, j, k)],
                    ids[(l + 1, j, k)],
                    ids[(l + 1, j + 1, k)],
                    ids[(l, j + 1, k)],
                    ids[(l, j, k + 1)],
                    ids[(l + 1, j, k + 1)],
                    ids[(l + 1, j + 1, k + 1)],
                    ids[(l, j + 1, k + 1)],
                ]
                add_elem(pid_capa(l), conn)


def bloque_y():
    ids = {}

    for l in range(capas_radiales + 1):
        for i, x in enumerate(xs):
            for k, z in enumerate(zs):
                ids[(l, i, k)] = add_nodo(punto_radial(x, a, z, l))

    for l in range(capas_radiales):
        for i in range(2 * d):
            for k in range(d):
                conn = [
                    ids[(l, i, k)],
                    ids[(l, i + 1, k)],
                    ids[(l + 1, i + 1, k)],
                    ids[(l + 1, i, k)],
                    ids[(l, i, k + 1)],
                    ids[(l, i + 1, k + 1)],
                    ids[(l + 1, i + 1, k + 1)],
                    ids[(l + 1, i, k + 1)],
                ]
                add_elem(pid_capa(l), conn)


def bloque_z():
    ids = {}

    for l in range(capas_radiales + 1):
        for i, x in enumerate(xs):
            for j, y in enumerate(ys):
                ids[(l, i, j)] = add_nodo(punto_radial(x, y, -a, l))

    for l in range(capas_radiales):
        for i in range(2 * d):
            for j in range(d):
                conn = [
                    ids[(l, i, j)],
                    ids[(l, i + 1, j)],
                    ids[(l, i + 1, j + 1)],
                    ids[(l, i, j + 1)],
                    ids[(l + 1, i, j)],
                    ids[(l + 1, i + 1, j)],
                    ids[(l + 1, i + 1, j + 1)],
                    ids[(l + 1, i, j + 1)],
                ]
                add_elem(pid_capa(l), conn)


bloque_x(-1.0)
bloque_x(1.0)
bloque_y()
bloque_z()


# Archivo LS-DYNA
salida.parent.mkdir(parents=True, exist_ok=True)
lineas = []

lineas.append("*KEYWORD\n")
lineas.extend(sec_solid())
lineas.extend(
    mat_block(
        mid_cito,
        mat["ro_citoplasma"],
        mat["c1_citoplasma"],
        mat["xk_citoplasma"],
        "Material_citoplasma",
    )
)
lineas.extend(
    mat_block(
        mid_nuc,
        mat["ro_nucleo"],
        mat["c1_nucleo"],
        mat["xk_nucleo"],
        "Material_nucleo",
    )
)
lineas.extend(part(pid_cito, mid_cito, "Citoplasma"))
lineas.extend(part(pid_nuc, mid_nuc, "Nucleo"))

lineas.append("$\n")
lineas.append("*NODE\n")

for nid in sorted(nodos):
    lineas.append(nodo_txt(nid, nodos[nid]))

lineas.append("$\n")
lineas.append("*ELEMENT_SOLID\n")

for eid, pid, conn in elems:
    lineas.append(elem_txt(eid, pid, conn))

lineas.append("*END\n")

with open(salida, "w", encoding="utf-8") as f:
    f.writelines(lineas)


# Resumen
n_nuc = sum(1 for eid, pid, conn in elems if pid == pid_nuc)
n_cito = sum(1 for eid, pid, conn in elems if pid == pid_cito)

print("Archivo:", salida)
print("Centro:", (cx, cy, cz))
print("Rc:", Rc)
print("eta:", eta)
print("Rn:", Rn)
print("q_cubo:", q_cubo)
print("Capas nucleo fuera del cubo:", capas_nuc)
print("Capas citoplasma radiales:", capas_cito)
print("Nodos:", len(nodos))
print("Elementos:", len(elems))
print("Citoplasma:", n_cito)
print("Nucleo:", n_nuc)
print("Material:", mat["preset"])
print("Ec kPa:", mat["Ec_kPa"])
print("En kPa:", mat["En_kPa"])
