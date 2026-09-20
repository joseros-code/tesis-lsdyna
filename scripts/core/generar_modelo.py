# generar_modelo.py
# ------------------------------------------------------------
# Genera solo el modelo completo.
# La celula se crea como archivo temporal y se borra al final.
#
# Ejemplo:
#   py scripts/core/generar_modelo.py --d 18 --r 0.005 --eta 0.67 --cz 0.01
# ------------------------------------------------------------

import argparse
import subprocess
import sys
from pathlib import Path

from definir_materiales import args_mat, leer_mat


core = Path(__file__).resolve().parent
base = Path(__file__).resolve().parents[2]


def leer_args():
    p = argparse.ArgumentParser(description="Genera el modelo completo.")
    p.add_argument("--densidad", "--d", dest="densidad", type=int, default=12)
    p.add_argument("--Rc", "--r", dest="Rc", type=float, default=0.005)
    p.add_argument("--eta", type=float, default=0.67)
    p.add_argument("--cx", type=float, default=0.015)
    p.add_argument("--cy", type=float, default=0.0)
    p.add_argument("--cz", type=float, default=0.01)
    args_mat(p)
    p.add_argument(
        "--modelo-base",
        "--base",
        dest="modelo_base",
        default="",
        help="Modelo base de Carlos. Si se omite usa originales/MallaCelula_DD.k.",
    )
    p.add_argument(
        "--prefijo-salida",
        "--tag-out",
        dest="prefijo_salida",
        default="",
        help="Prefijo opcional para los archivos generados.",
    )

    return p.parse_args()


def tag(valor, decimales):
    return f"{valor:.{decimales}f}".replace(".", "p")


def buscar_base(densidad):
    nombre = f"MallaCelula_{densidad:02d}.k"
    ruta_directa = base / "originales" / nombre

    if ruta_directa.exists():
        return ruta_directa

    candidatos = sorted((base / "originales").rglob(nombre))

    if candidatos:
        return candidatos[0]

    raise FileNotFoundError(f"No se encontro {nombre} en {base / 'originales'}")


def run(cmd):
    print("", flush=True)
    print("Ejecutando:", flush=True)
    print(" ".join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True)


args = leer_args()
mat = leer_mat(args)

d_tag = f"D{args.densidad:02d}"
r_tag = "Rc" + tag(args.Rc, 4)
eta_tag = "eta" + tag(args.eta, 3)
rel_tag = "EnEc" + tag(mat["relacion_En_Ec"], 2)
mat_tag = args.material_preset

if args.prefijo_salida:
    prefijo = args.prefijo_salida
else:
    prefijo = f"{d_tag}_{r_tag}_{eta_tag}_{mat_tag}_{rel_tag}"

if args.modelo_base:
    modelo_base = Path(args.modelo_base)

    if not modelo_base.is_absolute():
        modelo_base = base / modelo_base
else:
    modelo_base = buscar_base(args.densidad)

archivo_celula = base / "generados" / f"_tmp_cell_{prefijo}.k"
archivo_modelo = base / "generados" / f"Modelo_{prefijo}.k"

if not modelo_base.exists():
    raise FileNotFoundError(f"No existe el modelo base: {modelo_base}")

run(
    [
        sys.executable,
        str(core / "generar_malla.py"),
        "--r",
        str(args.Rc),
        "--eta",
        str(args.eta),
        "--d",
        str(args.densidad),
        "--cx",
        str(args.cx),
        "--cy",
        str(args.cy),
        "--cz",
        str(args.cz),
        "--mat",
        args.material_preset,
        "--Ec",
        str(mat["Ec_kPa"]),
        "--En",
        str(mat["En_kPa"]),
        "--nu-c",
        str(args.nu_citoplasma),
        "--nu-n",
        str(args.nu_nucleo),
        "--ro-c",
        str(args.ro_citoplasma),
        "--ro-n",
        str(args.ro_nucleo),
        "--out",
        str(archivo_celula),
    ]
)

run(
    [
        sys.executable,
        str(core / "integrar_celula.py"),
        "--mod",
        str(modelo_base),
        "--cell",
        str(archivo_celula),
        "--out",
        str(archivo_modelo),
        "--mat",
        args.material_preset,
        "--Ec",
        str(mat["Ec_kPa"]),
        "--En",
        str(mat["En_kPa"]),
        "--nu-c",
        str(args.nu_citoplasma),
        "--nu-n",
        str(args.nu_nucleo),
        "--ro-c",
        str(args.ro_citoplasma),
        "--ro-n",
        str(args.ro_nucleo),
    ]
)


if archivo_celula.exists():
    archivo_celula.unlink()

print("terminado.")
print("Modelo:", archivo_modelo)
