# barrido_gamma.py
# ------------------------------------------------------------
# Barrido de contraste mecanico nucleo-citoplasma para LS-DYNA.
#
# Geometria fija de la campania:
#   D18, Rc=0.005 mm, eta=0.67, centro=(0.015, 0.0, 0.01) mm
#
# La variable del barrido es:
#   Gamma = En / Ec
#   Ec = En / Gamma
#
# Ejemplos de uso:
#   py scripts/campaigns/barrido_gamma.py --dry-run
#   py scripts/campaigns/barrido_gamma.py --max-priority 3
#   py scripts/campaigns/barrido_gamma.py --only NC-07 NC-08 --dt2ms 0
#   py scripts/campaigns/barrido_gamma.py
# ------------------------------------------------------------

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


scripts_root = Path(__file__).resolve().parent.parent
core = scripts_root / "core"
base = scripts_root.parent


D = 18
Rc = 0.005
eta = 0.67
cx = 0.015
cy = 0.0
cz = 0.01
nu_c = 0.49999
nu_n = 0.49999
rho_c = 1.4e-9
rho_n = 1.4e-9


casos = [
    {
        "prioridad": 1,
        "case_id": "NC-07",
        "En_kPa": 300.0,
        "Gamma": 1.5,
        "Ec_tabla": 200.000000,
    },
    {
        "prioridad": 2,
        "case_id": "NC-08",
        "En_kPa": 300.0,
        "Gamma": 2.0,
        "Ec_tabla": 150.000000,
    },
    {
        "prioridad": 3,
        "case_id": "NC-09",
        "En_kPa": 300.0,
        "Gamma": 3.0,
        "Ec_tabla": 100.000000,
    },
    {
        "prioridad": 4,
        "case_id": "NC-04",
        "En_kPa": 225.0,
        "Gamma": 1.5,
        "Ec_tabla": 150.000000,
    },
    {
        "prioridad": 5,
        "case_id": "NC-05",
        "En_kPa": 225.0,
        "Gamma": 2.0,
        "Ec_tabla": 112.500000,
    },
    {
        "prioridad": 6,
        "case_id": "NC-06",
        "En_kPa": 225.0,
        "Gamma": 3.0,
        "Ec_tabla": 75.000000,
    },
    {
        "prioridad": 7,
        "case_id": "NC-01",
        "En_kPa": 100.0,
        "Gamma": 1.5,
        "Ec_tabla": 66.666667,
    },
    {
        "prioridad": 8,
        "case_id": "NC-02",
        "En_kPa": 100.0,
        "Gamma": 2.0,
        "Ec_tabla": 50.000000,
    },
    {
        "prioridad": 9,
        "case_id": "NC-03",
        "En_kPa": 100.0,
        "Gamma": 3.0,
        "Ec_tabla": 33.333333,
    },
]


def leer_args():
    p = argparse.ArgumentParser(
        description="Genera barrido de contraste mecanico D18."
    )
    p.add_argument(
        "--only",
        nargs="+",
        default=[],
        help="Genera solo los casos indicados, por ejemplo: --only NC-07 NC-08.",
    )
    p.add_argument(
        "--max-priority",
        type=int,
        default=None,
        help="Genera casos con prioridad menor o igual al valor indicado.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Imprime comandos y archivos de trazabilidad, sin generar modelos.",
    )
    p.add_argument(
        "--out-dir",
        default="generados/barrido_contraste_D18",
        help="Carpeta de salida.",
    )
    p.add_argument(
        "--dt2ms",
        type=float,
        default=-1.11e-10,
        help="DT2MS para escalado automatico de masa. Use 0 para apagarlo.",
    )
    p.add_argument(
        "--modelo-base",
        "--base",
        dest="modelo_base",
        default="",
        help="Modelo base D18. Si se omite busca MallaCelula_18.k.",
    )

    return p.parse_args()


def tag_num(valor, decimales=6):
    if abs(valor - round(valor)) < 1e-9:
        texto = str(int(round(valor)))
    else:
        texto = f"{valor:.{decimales}f}".rstrip("0").rstrip(".")

    return texto.replace(".", "p")


def nombre_modelo(caso):
    case_id = caso["case_id"]
    En = tag_num(caso["En_kPa"])
    Ec = tag_num(caso["Ec_kPa"])
    Gamma = f"{caso['Gamma']:.2f}".replace(".", "p")
    return f"modelo_{case_id}_En{En}_Ec{Ec}_Gamma{Gamma}.k"


def buscar_base():
    originales = base / "originales"

    if not originales.exists():
        raise FileNotFoundError(f"No existe la carpeta de originales: {originales}")

    nombre = "MallaCelula_18.k"
    directa = originales / nombre

    if directa.exists():
        return directa

    candidatos = sorted(originales.rglob(nombre))

    if candidatos:
        return candidatos[0]

    raise FileNotFoundError(f"No se encontro {nombre} dentro de {originales}")


def preparar_modelo_base(ruta_arg):
    if ruta_arg:
        ruta = Path(ruta_arg)

        if not ruta.is_absolute():
            ruta = base / ruta
    else:
        ruta = buscar_base()

    if not ruta.parent.exists():
        raise FileNotFoundError(f"No existe la carpeta del modelo base: {ruta.parent}")

    if not ruta.exists():
        raise FileNotFoundError(f"No existe el modelo base: {ruta}")

    return ruta


def validar_casos(lista):
    ids = set()
    tol = 1e-6

    for caso in lista:
        case_id = caso["case_id"]

        if case_id in ids:
            raise ValueError(f"Identificador de caso duplicado: {case_id}")

        ids.add(case_id)

        En = caso["En_kPa"]
        Gamma = caso["Gamma"]
        Ec = En / Gamma

        if En <= 0.0:
            raise ValueError(f"{case_id}: En_kPa debe ser positivo")

        if Gamma <= 1.0:
            raise ValueError(f"{case_id}: Gamma debe ser mayor que 1")

        if Ec <= 0.0:
            raise ValueError(f"{case_id}: Ec_kPa debe ser positivo")

        if not math.isclose(Ec, caso["Ec_tabla"], rel_tol=0.0, abs_tol=tol):
            raise ValueError(
                f"{case_id}: Ec calculado no coincide con la tabla "
                f"({Ec:.9g} != {caso['Ec_tabla']:.9g})"
            )

        caso["Ec_kPa"] = Ec
        caso["archivo_modelo"] = nombre_modelo(caso)


def seleccionar_casos(args):
    seleccion = {caso["case_id"] for caso in casos}

    if args.only:
        pedidos = {valor.upper() for valor in args.only}
        conocidos = {caso["case_id"] for caso in casos}
        desconocidos = sorted(pedidos - conocidos)

        if desconocidos:
            raise ValueError("Casos no reconocidos: " + ", ".join(desconocidos))

        seleccion = pedidos

    if args.max_priority is not None:
        if args.max_priority < 1:
            raise ValueError("--max-priority debe ser mayor o igual a 1")

        por_prioridad = {
            caso["case_id"]
            for caso in casos
            if caso["prioridad"] <= args.max_priority
        }
        seleccion = seleccion & por_prioridad

    for caso in casos:
        if caso["case_id"] in seleccion:
            caso["estado"] = "dry-run" if args.dry_run else "pendiente"
        else:
            caso["estado"] = "omitido_por_filtro"

    return [caso for caso in casos if caso["case_id"] in seleccion]


def run(cmd, dry_run):
    print("")
    print(" ".join(str(c) for c in cmd), flush=True)

    if not dry_run:
        subprocess.run(cmd, check=True)


def escribir_csv(out_dir, dt2ms):
    ruta = out_dir / "casos_contraste_material.csv"
    campos = [
        "prioridad",
        "case_id",
        "En_kPa",
        "Gamma_En_Ec",
        "Ec_kPa",
        "D",
        "Rc_mm",
        "eta",
        "cx_mm",
        "cy_mm",
        "cz_mm",
        "nu_c",
        "nu_n",
        "rho_c_ton_mm3",
        "rho_n_ton_mm3",
        "dt2ms",
        "archivo_modelo",
        "estado",
    ]

    with open(ruta, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=campos)
        writer.writeheader()

        for caso in casos:
            writer.writerow(
                {
                    "prioridad": caso["prioridad"],
                    "case_id": caso["case_id"],
                    "En_kPa": f"{caso['En_kPa']:.6f}",
                    "Gamma_En_Ec": f"{caso['Gamma']:.6f}",
                    "Ec_kPa": f"{caso['Ec_kPa']:.6f}",
                    "D": D,
                    "Rc_mm": f"{Rc:.6f}",
                    "eta": f"{eta:.6f}",
                    "cx_mm": f"{cx:.6f}",
                    "cy_mm": f"{cy:.6f}",
                    "cz_mm": f"{cz:.6f}",
                    "nu_c": f"{nu_c:.5f}",
                    "nu_n": f"{nu_n:.5f}",
                    "rho_c_ton_mm3": f"{rho_c:.6E}",
                    "rho_n_ton_mm3": f"{rho_n:.6E}",
                    "dt2ms": f"{dt2ms:.6E}",
                    "archivo_modelo": caso["archivo_modelo"],
                    "estado": caso["estado"],
                }
            )

    return ruta


def tabla_markdown(lista):
    lineas = [
        "| Prioridad | Case ID | En [kPa] | Gamma = En/Ec | Ec [kPa] | Estado |",
        "| --------: | ------- | -------: | -------------: | -------: | ------ |",
    ]

    for caso in lista:
        lineas.append(
            "| {prioridad} | {case_id} | {En:.6f} | {Gamma:.6f} | "
            "{Ec:.6f} | {estado} |".format(
                prioridad=caso["prioridad"],
                case_id=caso["case_id"],
                En=caso["En_kPa"],
                Gamma=caso["Gamma"],
                Ec=caso["Ec_kPa"],
                estado=caso["estado"],
            )
        )

    return "\n".join(lineas)


def escribir_txt(out_dir):
    ruta = out_dir / "casos_contraste_material.txt"
    texto = [
        "# Barrido de contraste material D18",
        "",
        "Este barrido evalua el contraste mecanico nucleo--citoplasma "
        "Gamma = En/Ec",
        "para tres niveles de rigidez nuclear. No incluye casos homogeneos.",
        "",
        tabla_markdown(casos),
        "",
    ]

    ruta.write_text("\n".join(texto), encoding="utf-8")
    return ruta


def escribir_orden(out_dir):
    ruta = out_dir / "orden_simulacion.txt"
    lineas = [
        "Orden recomendado de simulacion",
        "",
        "Simular primero los casos de prioridad 1 a 3.",
        "No interpretar la generacion exitosa del archivo .k como validacion fisica.",
        "Revisar d3hsp, glstat, matsum, nodout y d3plot despues de cada simulacion.",
        "",
    ]

    for caso in casos:
        lineas.append(
            "{prioridad}. {case_id}: En={En:.6f} kPa, Ec={Ec:.6f} kPa, "
            "Gamma={Gamma:.6f}, estado={estado}".format(
                prioridad=caso["prioridad"],
                case_id=caso["case_id"],
                En=caso["En_kPa"],
                Ec=caso["Ec_kPa"],
                Gamma=caso["Gamma"],
                estado=caso["estado"],
            )
        )

    ruta.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    return ruta


def escribir_trazabilidad(out_dir, dt2ms):
    out_dir.mkdir(parents=True, exist_ok=True)
    rutas = [
        escribir_csv(out_dir, dt2ms),
        escribir_txt(out_dir),
        escribir_orden(out_dir),
    ]

    return rutas


def imprimir_tabla(lista):
    print("")
    print("Casos de la campania:")
    print("Prioridad  Caso    En[kPa]     Gamma      Ec[kPa]      Estado")
    print("---------  ------  ----------  ---------  -----------  ----------------")

    for caso in lista:
        print(
            f"{caso['prioridad']:>9}  "
            f"{caso['case_id']:<6}  "
            f"{caso['En_kPa']:>10.6f}  "
            f"{caso['Gamma']:>9.6f}  "
            f"{caso['Ec_kPa']:>11.6f}  "
            f"{caso['estado']}"
        )


def main():
    args = leer_args()
    validar_casos(casos)

    modelo_base = preparar_modelo_base(args.modelo_base)
    out_dir = Path(args.out_dir)

    if not out_dir.is_absolute():
        out_dir = base / out_dir

    seleccionados = seleccionar_casos(args)

    if not seleccionados:
        raise ValueError("No hay casos seleccionados para generar")

    print("Barrido de contraste material")
    print("Modelo base:", modelo_base)
    print("Salida:", out_dir)
    print("Densidad:", D)
    print("Centro:", (cx, cy, cz))
    print("Rc:", Rc)
    print("eta:", eta)
    print("DT2MS:", args.dt2ms)
    print("Casos seleccionados:", len(seleccionados))

    out_dir.mkdir(parents=True, exist_ok=True)

    for caso in seleccionados:
        case_id = caso["case_id"]
        Ec_kPa = caso["Ec_kPa"]
        En_kPa = caso["En_kPa"]
        tmp_cell = out_dir / f"_tmp_cell_{case_id}.k"
        modelo_out = out_dir / caso["archivo_modelo"]

        print("")
        print("===", case_id, "===")
        print("Prioridad:", caso["prioridad"])
        print("En kPa:", f"{En_kPa:.6f}")
        print("Gamma En/Ec:", f"{caso['Gamma']:.6f}")
        print("Ec kPa:", f"{Ec_kPa:.6f}")

        cmd_mesh = [
            sys.executable,
            str(core / "generar_malla.py"),
            "--d",
            str(D),
            "--r",
            str(Rc),
            "--eta",
            str(eta),
            "--cx",
            str(cx),
            "--cy",
            str(cy),
            "--cz",
            str(cz),
            "--mat",
            "manual",
            "--Ec",
            str(Ec_kPa),
            "--En",
            str(En_kPa),
            "--nu-c",
            str(nu_c),
            "--nu-n",
            str(nu_n),
            "--ro-c",
            str(rho_c),
            "--ro-n",
            str(rho_n),
            "--out",
            str(tmp_cell),
        ]

        cmd_add = [
            sys.executable,
            str(core / "integrar_celula.py"),
            "--mod",
            str(modelo_base),
            "--cell",
            str(tmp_cell),
            "--out",
            str(modelo_out),
            "--mat",
            "manual",
            "--Ec",
            str(Ec_kPa),
            "--En",
            str(En_kPa),
            "--nu-c",
            str(nu_c),
            "--nu-n",
            str(nu_n),
            "--ro-c",
            str(rho_c),
            "--ro-n",
            str(rho_n),
            "--dt2ms",
            str(args.dt2ms),
        ]

        run(cmd_mesh, args.dry_run)
        run(cmd_add, args.dry_run)

        if args.dry_run:
            print("Borrar temporal:", tmp_cell)
        elif tmp_cell.exists():
            tmp_cell.unlink()

        if not args.dry_run:
            caso["estado"] = "generado"

    rutas = escribir_trazabilidad(out_dir, args.dt2ms)

    imprimir_tabla(casos)

    print("")
    print("Archivos de trazabilidad:")

    for ruta in rutas:
        print(" ", ruta)

    print("")
    print("Ejemplos de ejecucion:")
    print("  py scripts/campaigns/barrido_gamma.py --dry-run")
    print("  py scripts/campaigns/barrido_gamma.py --max-priority 3")
    print(
        "  py scripts/campaigns/barrido_gamma.py "
        "--only NC-07 NC-08 --dt2ms 0"
    )
    print("  py scripts/campaigns/barrido_gamma.py")

    print("")
    print("Barrido de contraste material terminado.")


if __name__ == "__main__":
    main()
