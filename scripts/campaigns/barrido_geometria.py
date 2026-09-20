# barrido_geometria.py
# ------------------------------------------------------------
# Barrido geometrico para la celula multicomponente en LS-DYNA.
#
# Mantiene fija la malla D18, el centro de la celula y la materialidad
# base NC-08:
#   Ec = 150 kPa, En = 300 kPa, Gamma = 2
#
# La campania varia:
#   Rc  = radio celular externo
#   eta = Rn / Rc
#
# Ejemplos:
#   py scripts/campaigns/barrido_geometria.py --dry-run
#   py scripts/campaigns/barrido_geometria.py --only GEO-01 GEO-05 GEO-09
#   py scripts/campaigns/barrido_geometria.py --dt2ms 0
#   py scripts/campaigns/barrido_geometria.py
# ------------------------------------------------------------

import argparse
import csv
import subprocess
import sys
from pathlib import Path


scripts_root = Path(__file__).resolve().parent.parent
core = scripts_root / "core"
base = scripts_root.parent


D = 18
cx = 0.015
cy = 0.0
cz = 0.01

Ec_kPa = 150.0
En_kPa = 300.0
Gamma = En_kPa / Ec_kPa
nu_c = 0.49999
nu_n = 0.49999
rho_c = 1.4e-9
rho_n = 1.4e-9


casos = [
    {"prioridad": 1, "case_id": "GEO-01", "Rc": 0.0045, "eta": 0.55},
    {"prioridad": 2, "case_id": "GEO-02", "Rc": 0.0045, "eta": 0.67},
    {"prioridad": 3, "case_id": "GEO-03", "Rc": 0.0045, "eta": 0.80},
    {"prioridad": 4, "case_id": "GEO-04", "Rc": 0.0050, "eta": 0.55},
    {"prioridad": 5, "case_id": "GEO-05", "Rc": 0.0050, "eta": 0.67},
    {"prioridad": 6, "case_id": "GEO-06", "Rc": 0.0050, "eta": 0.80},
    {"prioridad": 7, "case_id": "GEO-07", "Rc": 0.0055, "eta": 0.55},
    {"prioridad": 8, "case_id": "GEO-08", "Rc": 0.0055, "eta": 0.67},
    {"prioridad": 9, "case_id": "GEO-09", "Rc": 0.0055, "eta": 0.80},
]


def leer_args():
    p = argparse.ArgumentParser(description="Genera barrido geometrico D18.")
    p.add_argument(
        "--only",
        nargs="+",
        default=[],
        help="Genera solo los casos indicados, por ejemplo: --only GEO-01 GEO-05.",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Imprime comandos sin generar modelos.",
    )
    p.add_argument(
        "--out-dir",
        default="generados/barrido_geometrico_D18",
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
    texto = f"{valor:.{decimales}f}".rstrip("0").rstrip(".")
    return texto.replace(".", "p")


def nombre_modelo(caso):
    return (
        f"modelo_{caso['case_id']}_"
        f"Rc{tag_num(caso['Rc'], 4)}_"
        f"eta{tag_num(caso['eta'], 2)}.k"
    )


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


def validar_casos():
    ids = set()

    for caso in casos:
        case_id = caso["case_id"]

        if case_id in ids:
            raise ValueError(f"Identificador de caso duplicado: {case_id}")

        ids.add(case_id)

        if caso["Rc"] <= 0.0:
            raise ValueError(f"{case_id}: Rc debe ser positivo")

        if caso["eta"] <= 0.0 or caso["eta"] >= 1.0:
            raise ValueError(f"{case_id}: eta debe estar entre 0 y 1")

        caso["Rn"] = caso["Rc"] * caso["eta"]
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
    ruta = out_dir / "casos_geometricos.csv"
    campos = [
        "prioridad",
        "case_id",
        "Rc_mm",
        "eta",
        "Rn_mm",
        "D",
        "cx_mm",
        "cy_mm",
        "cz_mm",
        "Ec_kPa",
        "En_kPa",
        "Gamma_En_Ec",
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
                    "Rc_mm": f"{caso['Rc']:.6f}",
                    "eta": f"{caso['eta']:.6f}",
                    "Rn_mm": f"{caso['Rn']:.6f}",
                    "D": D,
                    "cx_mm": f"{cx:.6f}",
                    "cy_mm": f"{cy:.6f}",
                    "cz_mm": f"{cz:.6f}",
                    "Ec_kPa": f"{Ec_kPa:.6f}",
                    "En_kPa": f"{En_kPa:.6f}",
                    "Gamma_En_Ec": f"{Gamma:.6f}",
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


def tabla_markdown():
    lineas = [
        "| Prioridad | Caso | Rc [mm] | eta | Rn [mm] | Estado |",
        "| --------: | ---- | ------: | --: | ------: | ------ |",
    ]

    for caso in casos:
        lineas.append(
            "| {prioridad} | {case_id} | {Rc:.6f} | {eta:.6f} | "
            "{Rn:.6f} | {estado} |".format(
                prioridad=caso["prioridad"],
                case_id=caso["case_id"],
                Rc=caso["Rc"],
                eta=caso["eta"],
                Rn=caso["Rn"],
                estado=caso["estado"],
            )
        )

    return "\n".join(lineas)


def escribir_txt(out_dir):
    ruta = out_dir / "casos_geometricos.txt"
    texto = [
        "# Barrido geometrico D18",
        "",
        "Este barrido evalua el efecto del radio celular externo Rc y de la "
        "relacion geometrica nuclear eta = Rn/Rc.",
        "",
        "La materialidad se mantiene fija como caso base NC-08: "
        "Ec=150 kPa, En=300 kPa, Gamma=2.",
        "",
        tabla_markdown(),
        "",
    ]

    ruta.write_text("\n".join(texto), encoding="utf-8")
    return ruta


def escribir_orden(out_dir):
    ruta = out_dir / "orden_geometrico.txt"
    lineas = [
        "Orden recomendado de simulacion geometrica",
        "",
        "Simular primero GEO-05 como caso base geometrico.",
        "Luego comparar variaciones de Rc con eta fijo y variaciones de eta con Rc fijo.",
        "No interpretar la generacion exitosa del archivo .k como validacion fisica.",
        "Revisar d3hsp, glstat, matsum, nodout y d3plot despues de cada simulacion.",
        "",
    ]

    for caso in casos:
        lineas.append(
            "{prioridad}. {case_id}: Rc={Rc:.6f} mm, eta={eta:.6f}, "
            "Rn={Rn:.6f} mm, estado={estado}".format(
                prioridad=caso["prioridad"],
                case_id=caso["case_id"],
                Rc=caso["Rc"],
                eta=caso["eta"],
                Rn=caso["Rn"],
                estado=caso["estado"],
            )
        )

    ruta.write_text("\n".join(lineas) + "\n", encoding="utf-8")
    return ruta


def escribir_trazabilidad(out_dir, dt2ms):
    out_dir.mkdir(parents=True, exist_ok=True)
    return [
        escribir_csv(out_dir, dt2ms),
        escribir_txt(out_dir),
        escribir_orden(out_dir),
    ]


def imprimir_tabla():
    print("")
    print("Casos geometricos:")
    print("Prioridad  Caso    Rc[mm]    eta       Rn[mm]    Estado")
    print("---------  ------  --------  --------  --------  ----------------")

    for caso in casos:
        print(
            f"{caso['prioridad']:>9}  "
            f"{caso['case_id']:<6}  "
            f"{caso['Rc']:>8.6f}  "
            f"{caso['eta']:>8.6f}  "
            f"{caso['Rn']:>8.6f}  "
            f"{caso['estado']}"
        )


def main():
    args = leer_args()
    validar_casos()

    modelo_base = preparar_modelo_base(args.modelo_base)
    out_dir = Path(args.out_dir)

    if not out_dir.is_absolute():
        out_dir = base / out_dir

    seleccionados = seleccionar_casos(args)

    if not seleccionados:
        raise ValueError("No hay casos seleccionados para generar")

    print("Barrido geometrico")
    print("Modelo base:", modelo_base)
    print("Salida:", out_dir)
    print("Densidad:", D)
    print("Centro:", (cx, cy, cz))
    print("Material base: NC-08")
    print("Ec kPa:", Ec_kPa)
    print("En kPa:", En_kPa)
    print("Gamma En/Ec:", Gamma)
    print("DT2MS:", args.dt2ms)
    print("Casos seleccionados:", len(seleccionados))

    out_dir.mkdir(parents=True, exist_ok=True)

    for caso in seleccionados:
        case_id = caso["case_id"]
        tmp_cell = out_dir / f"_tmp_cell_{case_id}.k"
        modelo_out = out_dir / caso["archivo_modelo"]

        print("")
        print("===", case_id, "===")
        print("Prioridad:", caso["prioridad"])
        print("Rc mm:", f"{caso['Rc']:.6f}")
        print("eta:", f"{caso['eta']:.6f}")
        print("Rn mm:", f"{caso['Rn']:.6f}")

        cmd_mesh = [
            sys.executable,
            str(core / "generar_malla.py"),
            "--d",
            str(D),
            "--r",
            str(caso["Rc"]),
            "--eta",
            str(caso["eta"]),
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

    imprimir_tabla()

    print("")
    print("Archivos de trazabilidad:")

    for ruta in rutas:
        print(" ", ruta)

    print("")
    print("Barrido geometrico terminado.")


if __name__ == "__main__":
    main()
