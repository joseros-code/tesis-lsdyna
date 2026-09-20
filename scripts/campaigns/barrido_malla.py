# barrido_malla.py
# ------------------------------------------------------------
# Genera varios modelos para sensibilidad de malla.
#
# Ejemplo:
#   py scripts/campaigns/barrido_malla.py --ds 8 10 12 14 16 18 20 22 24 --dt2ms -1.110e-10
#
# Cada densidad usa por defecto el modelo base:
#   originales/MallaCelula_DD.k
#
# y genera solo:
#   generados/Modelo_sens_malla_DDD.k
# ------------------------------------------------------------

import argparse
import subprocess
import sys
from pathlib import Path


scripts_root = Path(__file__).resolve().parent.parent
core = scripts_root / "core"
base = scripts_root.parent

if str(core) not in sys.path:
    sys.path.insert(0, str(core))

from definir_materiales import args_mat, leer_mat


def leer_args():
    p = argparse.ArgumentParser(description="Genera una familia de mallas.")
    p.add_argument(
        "--densidades",
        "--ds",
        dest="densidades",
        type=int,
        nargs="+",
        default=[8, 10, 12, 14, 16, 18, 20, 22, 24],
        help="Lista de densidades de la malla celular.",
    )
    p.add_argument("--Rc", "--r", dest="Rc", type=float, default=0.005)
    p.add_argument("--eta", type=float, default=0.67)
    p.add_argument("--cx", type=float, default=0.015)
    p.add_argument("--cy", type=float, default=0.0)
    p.add_argument("--cz", type=float, default=0.01)
    p.add_argument(
        "--dt2ms",
        type=float,
        default=-1.110e-10,
        help="DT2MS para escalado automatico de masa. Use 0 para apagarlo.",
    )
    p.add_argument(
        "--tag",
        default="sens_malla",
        help="Etiqueta corta para agrupar los archivos generados.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Muestra comandos sin ejecutarlos.",
    )
    args_mat(p)

    return p.parse_args()


def buscar_base(densidad):
    nombre = f"MallaCelula_{densidad:02d}.k"
    ruta_directa = base / "originales" / nombre

    if ruta_directa.exists():
        return ruta_directa

    candidatos = sorted((base / "originales").rglob(nombre))

    if candidatos:
        return candidatos[0]

    return None


def run(cmd, dry_run):
    print("", flush=True)
    print(" ".join(str(c) for c in cmd), flush=True)

    if not dry_run:
        subprocess.run(cmd, check=True)


args = leer_args()
mat = leer_mat(args)

print("Lote de sensibilidad de malla")
print("Densidades:", args.densidades)
print("Rc:", args.Rc)
print("eta:", args.eta)
print("centro:", (args.cx, args.cy, args.cz))
print("DT2MS:", args.dt2ms)
print("material preset:", args.material_preset)
print("Ec kPa:", mat["Ec_kPa"])
print("En kPa:", mat["En_kPa"])
print("En/Ec:", mat["relacion_En_Ec"])

for densidad in args.densidades:
    modelo_base = buscar_base(densidad)

    if modelo_base is None:
        print("")
        print(
            "Saltando densidad",
            densidad,
            "porque no se encontro MallaCelula_"
            f"{densidad:02d}.k dentro de:",
            base / "originales",
        )
        continue

    celula_tmp = base / "generados" / f"_tmp_cell_D{densidad:02d}.k"
    modelo_out = base / "generados" / f"Modelo_{args.tag}_D{densidad:02d}.k"

    cmd_mesh = [
        sys.executable,
        str(core / "generar_malla.py"),
        "--d",
        str(densidad),
        "--r",
        str(args.Rc),
        "--eta",
        str(args.eta),
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
        str(celula_tmp),
    ]

    cmd_add = [
        sys.executable,
        str(core / "integrar_celula.py"),
        "--mod",
        str(modelo_base),
        "--cell",
        str(celula_tmp),
        "--out",
        str(modelo_out),
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

    if args.dt2ms != 0:
        cmd_add.extend(["--dt2ms", str(args.dt2ms)])

    run(cmd_mesh, args.dry_run)
    run(cmd_add, args.dry_run)

    if args.dry_run:
        print("Borrar temporal:", celula_tmp)
    elif celula_tmp.exists():
        celula_tmp.unlink()

print("")
print("Lote terminado.")
