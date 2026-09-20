#!/usr/bin/env python3
"""Extrae y resume las campañas paramétricas de la tesis.

El script usa LS-PrePost para leer los d3plot binarios, identifica de forma
reproducible los nodos frontales del citoplasma (PID 2) y del núcleo (PID 5),
y genera tablas y figuras normalizadas.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


TESIS_LSDYNA = Path(__file__).resolve().parents[2]
RESULTADOS = TESIS_LSDYNA / "resultados"
GENERADOS = TESIS_LSDYNA / "generados"
SALIDA = TESIS_LSDYNA / "resultados_procesados"
LSPP_PREDETERMINADO = Path(
    r"C:\Program Files\LS-DYNA Suite R16.1 Student\lspp\lsprepost4.12.exe"
)


@dataclass(frozen=True)
class Caso:
    campana: str
    codigo: str
    carpeta: str
    en_kpa: float | None = None
    ec_kpa: float | None = None
    gamma: float | None = None
    rc_mm: float | None = None
    eta: float | None = None


CASOS_MATERIAL = [
    Caso("material", "NC-01", "En100_Ec66p67_Gamma1p5", 100.0, 66.666667, 1.5, 0.0050, 0.67),
    Caso("material", "NC-02", "En100_Ec50_Gamma2", 100.0, 50.0, 2.0, 0.0050, 0.67),
    Caso("material", "NC-03", "En100_Ec33p3_Gamma3", 100.0, 33.333333, 3.0, 0.0050, 0.67),
    Caso("material", "NC-04", "En225_Ec150_Gamma1p5", 225.0, 150.0, 1.5, 0.0050, 0.67),
    Caso("material", "NC-05", "En225_Ec112p5_Gamma2p0", 225.0, 112.5, 2.0, 0.0050, 0.67),
    Caso("material", "NC-06", "En225_Ec75_Gamma3p0", 225.0, 75.0, 3.0, 0.0050, 0.67),
    Caso("material", "NC-07", "En300_Ec200_Gamma1p5", 300.0, 200.0, 1.5, 0.0050, 0.67),
    Caso("material", "NC-08", "En300_Ec150_Gamma2p0", 300.0, 150.0, 2.0, 0.0050, 0.67),
    Caso("material", "NC-09", "En300_Ec100_Gamma3p0", 300.0, 100.0, 3.0, 0.0050, 0.67),
]

CASOS_GEOMETRIA = [
    Caso("geometria", "GEO-01", "geo_01", rc_mm=0.0045, eta=0.55),
    Caso("geometria", "GEO-02", "geo_02", rc_mm=0.0045, eta=0.67),
    Caso("geometria", "GEO-03", "geo_03", rc_mm=0.0045, eta=0.80),
    Caso("geometria", "GEO-04", "geo_04", rc_mm=0.0050, eta=0.55),
    Caso("geometria", "GEO-05", "geo_05", rc_mm=0.0050, eta=0.67),
    Caso("geometria", "GEO-06", "geo_06", rc_mm=0.0050, eta=0.80),
    Caso("geometria", "GEO-07", "geo_07", rc_mm=0.0055, eta=0.55),
    Caso("geometria", "GEO-08", "geo_08", rc_mm=0.0055, eta=0.67),
    Caso("geometria", "GEO-09", "geo_09", rc_mm=0.0055, eta=0.80),
]


def carpeta_resultado(caso: Caso) -> Path:
    base = RESULTADOS / ("gamma" if caso.campana == "material" else "sens_geo")
    return base / caso.carpeta if caso.carpeta else base / "__ausente__"


def localizar_modelo(caso: Caso, resultado: Path) -> Path:
    modelos_locales = sorted(resultado.glob("*.k")) if resultado.is_dir() else []
    if modelos_locales:
        return modelos_locales[0]

    patron = f"modelo_{caso.codigo}*.k"
    coincidencias = sorted(GENERADOS.rglob(patron))
    if not coincidencias:
        raise FileNotFoundError(f"No se encontró el modelo de {caso.codigo}")
    return coincidencias[0]


def _campos_numericos(linea: str) -> list[str]:
    return linea.replace(",", " ").split()


def nodos_frontales(modelo: Path) -> dict[str, tuple[int, float]]:
    """Devuelve (nid, x0) para PID 2 y PID 5.

    El frente se define como el nodo perteneciente a elementos sólidos de la
    parte correspondiente que posee la mayor coordenada X inicial.
    """
    nodos_por_pid: dict[int, set[int]] = {2: set(), 5: set()}
    seccion = ""
    with modelo.open("r", encoding="latin-1", errors="ignore") as archivo:
        for linea in archivo:
            stripped = linea.strip()
            if not stripped or stripped.startswith("$"):
                continue
            if stripped.startswith("*"):
                seccion = stripped.upper()
                continue
            if seccion == "*ELEMENT_SOLID":
                partes = _campos_numericos(linea)
                if len(partes) < 3:
                    continue
                try:
                    pid = int(partes[1])
                except ValueError:
                    continue
                if pid in nodos_por_pid:
                    nodos_por_pid[pid].update(
                        int(valor) for valor in partes[2:] if int(valor) > 0
                    )

    mejores: dict[int, tuple[int, float]] = {}
    pendientes = nodos_por_pid[2] | nodos_por_pid[5]
    with modelo.open("r", encoding="latin-1", errors="ignore") as archivo:
        seccion = ""
        for linea in archivo:
            stripped = linea.strip()
            if not stripped or stripped.startswith("$"):
                continue
            if stripped.startswith("*"):
                seccion = stripped.upper()
                continue
            if seccion != "*NODE":
                continue
            partes = _campos_numericos(linea)
            if len(partes) < 4:
                continue
            try:
                nid = int(partes[0])
            except ValueError:
                continue
            if nid not in pendientes:
                continue
            x = float(partes[1].replace("D", "E"))
            for pid in (2, 5):
                if nid in nodos_por_pid[pid]:
                    anterior = mejores.get(pid)
                    if anterior is None or x > anterior[1]:
                        mejores[pid] = (nid, x)

    if 2 not in mejores or 5 not in mejores:
        raise RuntimeError(f"No se pudieron identificar ambos frentes en {modelo}")
    return {"citoplasma": mejores[2], "nucleo": mejores[5]}


def estado_simulacion(resultado: Path) -> str:
    d3hsp = resultado / "d3hsp"
    if not d3hsp.exists():
        return "sin_d3hsp"
    with d3hsp.open("rb") as archivo:
        tamano = d3hsp.stat().st_size
        archivo.seek(max(0, tamano - 500_000))
        lineas = archivo.read().decode("latin-1", errors="ignore").splitlines()
    for linea in reversed(lineas):
        compacta = re.sub(r"\s+", "", linea.lower())
        if compacta.startswith("normaltermination"):
            return "normal"
        if compacta.startswith("errortermination"):
            return "error"
    return "sin_marca_terminacion"


def escribir_cfile(
    ruta: Path,
    modelo: Path,
    nodos: dict[str, tuple[int, float]],
    salidas: dict[str, Path],
) -> None:
    lineas = [
        'open d3plot "d3plot"',
        "ac",
        f'open keyword "{modelo.name}"',
        "Message 1",
        "Message 0",
    ]
    for componente in ("citoplasma", "nucleo"):
        nid = nodos[componente][0]
        lineas.extend(
            [
                "genselect clear all",
                "genselect target node",
                "-M X1",
                f"genselect node add node {nid}/0",
                "ntime 1",
                f'xyplot 1 savefile ms_csv "{salidas[componente].as_posix()}" 1 all',
                "xyplot 1 donemenu",
                "deletewin 1",
            ]
        )
    lineas.append("quit")
    ruta.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def leer_csv_lspp(ruta: Path) -> list[tuple[float, float]]:
    datos: list[tuple[float, float]] = []
    with ruta.open("r", encoding="utf-8-sig", errors="ignore", newline="") as archivo:
        lector = csv.reader(archivo)
        next(lector, None)
        next(lector, None)
        for fila in lector:
            if len(fila) < 2 or not fila[0].strip() or not fila[1].strip():
                continue
            datos.append((float(fila[0]), float(fila[1])))
    if not datos:
        raise RuntimeError(f"LS-PrePost no generó datos válidos en {ruta}")
    return datos


def extraer_con_lspp(
    caso: Caso,
    resultado: Path,
    modelo: Path,
    nodos: dict[str, tuple[int, float]],
    lspp: Path,
    destino_crudo: Path,
) -> dict[str, list[tuple[float, float]]]:
    destino_crudo.mkdir(parents=True, exist_ok=True)
    salidas = {
        nombre: destino_crudo / f"{caso.codigo.lower()}_{nombre}_x.csv"
        for nombre in ("citoplasma", "nucleo")
    }
    if all(ruta.exists() for ruta in salidas.values()):
        return {nombre: leer_csv_lspp(ruta) for nombre, ruta in salidas.items()}

    with tempfile.TemporaryDirectory(prefix=f"{caso.codigo}_", dir=SALIDA) as tmp_nombre:
        tmp = Path(tmp_nombre)
        for d3plot in resultado.glob("d3plot*"):
            os.link(d3plot, tmp / d3plot.name)
        os.link(modelo, tmp / modelo.name)
        cfile = tmp / "extraer.cfile"
        escribir_cfile(cfile, modelo, nodos, salidas)
        proceso = subprocess.run(
            [str(lspp), f"c={cfile}"],
            cwd=tmp,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            timeout=180,
            check=False,
        )
        if proceso.returncode != 0:
            mensaje = proceso.stderr.strip() or f"código {proceso.returncode}"
            raise RuntimeError(f"LS-PrePost falló para {caso.codigo}: {mensaje}")
    return {nombre: leer_csv_lspp(ruta) for nombre, ruta in salidas.items()}


def interpolar(datos: list[tuple[float, float]], tiempo: float) -> float:
    if tiempo < datos[0][0] or tiempo > datos[-1][0]:
        return math.nan
    for (t0, x0), (t1, x1) in zip(datos, datos[1:]):
        if t0 <= tiempo <= t1:
            if t1 == t0:
                return x1
            fraccion = (tiempo - t0) / (t1 - t0)
            return x0 + fraccion * (x1 - x0)
    return datos[-1][1]


def escribir_tabla(ruta: Path, filas: Iterable[dict], campos: list[str]) -> None:
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with ruta.open("w", encoding="utf-8-sig", newline="") as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=campos)
        escritor.writeheader()
        escritor.writerows(filas)


def procesar_campana(
    campana: str,
    casos: list[Caso],
    lspp: Path,
    tiempo_ref_s: float,
) -> tuple[list[dict], list[dict]]:
    destino = SALIDA / campana
    crudos = destino / "series_lspp"
    resumen: list[dict] = []
    trayectorias: list[dict] = []

    for caso in casos:
        resultado = carpeta_resultado(caso)
        d3plot = resultado / "d3plot"
        if not d3plot.exists():
            print(f"[{caso.codigo}] pendiente: no existe d3plot", flush=True)
            resumen.append(
                {
                    "campana": campana,
                    "caso": caso.codigo,
                    "En_kPa": caso.en_kpa,
                    "Ec_kPa": caso.ec_kpa,
                    "Gamma": caso.gamma,
                    "Rc_mm": caso.rc_mm,
                    "eta": caso.eta,
                    "estado": "pendiente",
                    "nodo_frontal_citoplasma": "",
                    "nodo_frontal_nucleo": "",
                    "tiempo_final_ms": "",
                    "tiempo_referencia_ms": tiempo_ref_s * 1000,
                    "lambda_citoplasma_ref": "",
                    "lambda_nucleo_ref": "",
                    "rezago_nucleo_ref": "",
                }
            )
            continue

        modelo = localizar_modelo(caso, resultado)
        nodos = nodos_frontales(modelo)
        print(
            f"[{caso.codigo}] citoplasma={nodos['citoplasma'][0]}, "
            f"núcleo={nodos['nucleo'][0]}",
            flush=True,
        )
        series = extraer_con_lspp(caso, resultado, modelo, nodos, lspp, crudos)
        estado = estado_simulacion(resultado)
        tiempos_comunes = sorted(
            set(t for t, _ in series["citoplasma"])
            & set(t for t, _ in series["nucleo"])
        )
        nucleo_por_t = dict(series["nucleo"])
        citoplasma_por_t = dict(series["citoplasma"])
        rc = float(caso.rc_mm)
        x0_c = nodos["citoplasma"][1]
        x0_n = nodos["nucleo"][1]
        for tiempo in tiempos_comunes:
            lambda_c = (citoplasma_por_t[tiempo] - x0_c) / rc
            lambda_n = (nucleo_por_t[tiempo] - x0_n) / rc
            trayectorias.append(
                {
                    "campana": campana,
                    "caso": caso.codigo,
                    "En_kPa": caso.en_kpa,
                    "Ec_kPa": caso.ec_kpa,
                    "Gamma": caso.gamma,
                    "Rc_mm": caso.rc_mm,
                    "eta": caso.eta,
                    "estado": estado,
                    "tiempo_s": tiempo,
                    "tiempo_ms": tiempo * 1000,
                    "x_citoplasma_mm": citoplasma_por_t[tiempo],
                    "x_nucleo_mm": nucleo_por_t[tiempo],
                    "lambda_citoplasma": lambda_c,
                    "lambda_nucleo": lambda_n,
                    "rezago_nucleo": lambda_c - lambda_n,
                }
            )

        x_c_ref = interpolar(series["citoplasma"], tiempo_ref_s)
        x_n_ref = interpolar(series["nucleo"], tiempo_ref_s)
        lambda_c_ref = (x_c_ref - x0_c) / rc
        lambda_n_ref = (x_n_ref - x0_n) / rc
        resumen.append(
            {
                "campana": campana,
                "caso": caso.codigo,
                "En_kPa": caso.en_kpa,
                "Ec_kPa": caso.ec_kpa,
                "Gamma": caso.gamma,
                "Rc_mm": caso.rc_mm,
                "eta": caso.eta,
                "estado": estado,
                "nodo_frontal_citoplasma": nodos["citoplasma"][0],
                "nodo_frontal_nucleo": nodos["nucleo"][0],
                "tiempo_final_ms": min(
                    series["citoplasma"][-1][0], series["nucleo"][-1][0]
                )
                * 1000,
                "tiempo_referencia_ms": tiempo_ref_s * 1000,
                "lambda_citoplasma_ref": lambda_c_ref,
                "lambda_nucleo_ref": lambda_n_ref,
                "rezago_nucleo_ref": lambda_c_ref - lambda_n_ref,
            }
        )

    campos_resumen = [
        "campana",
        "caso",
        "En_kPa",
        "Ec_kPa",
        "Gamma",
        "Rc_mm",
        "eta",
        "estado",
        "nodo_frontal_citoplasma",
        "nodo_frontal_nucleo",
        "tiempo_final_ms",
        "tiempo_referencia_ms",
        "lambda_citoplasma_ref",
        "lambda_nucleo_ref",
        "rezago_nucleo_ref",
    ]
    campos_trayectorias = [
        "campana",
        "caso",
        "En_kPa",
        "Ec_kPa",
        "Gamma",
        "Rc_mm",
        "eta",
        "estado",
        "tiempo_s",
        "tiempo_ms",
        "x_citoplasma_mm",
        "x_nucleo_mm",
        "lambda_citoplasma",
        "lambda_nucleo",
        "rezago_nucleo",
    ]
    escribir_tabla(destino / f"resumen_{campana}.csv", resumen, campos_resumen)
    escribir_tabla(
        destino / f"trayectorias_{campana}.csv", trayectorias, campos_trayectorias
    )
    return resumen, trayectorias


def _numero(valor) -> float:
    if valor is None or valor == "":
        return math.nan
    return float(valor)


def generar_figuras(
    campana: str,
    resumen: list[dict],
    trayectorias: list[dict],
) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise RuntimeError(
            "Falta matplotlib. Instale las dependencias con "
            "`python -m pip install -r scripts/postprocess/requirements.txt`."
        ) from exc

    destino = SALIDA / campana / "figuras"
    destino.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "figure.dpi": 120,
            "savefig.dpi": 300,
        }
    )

    if campana == "material":
        paneles = [100.0, 225.0, 300.0]
        parametro_panel = "En_kPa"
        parametro_linea = "Gamma"
        colores = {1.5: "#0072B2", 2.0: "#E69F00", 3.0: "#009E73"}
        titulo_panel = lambda valor: rf"$E_n={valor:g}$ kPa"
        etiqueta_linea = lambda valor: rf"$\Gamma_E={valor:g}$"
        filas_heat = [100.0, 225.0, 300.0]
        columnas_heat = [1.5, 2.0, 3.0]
        fila_key, col_key = "En_kPa", "Gamma"
        etiqueta_y_heat = r"$E_n$ [kPa]"
        etiqueta_x_heat = r"$\Gamma_E=E_n/E_c$"
    else:
        paneles = [0.0045, 0.0050, 0.0055]
        parametro_panel = "Rc_mm"
        parametro_linea = "eta"
        colores = {0.55: "#0072B2", 0.67: "#E69F00", 0.80: "#009E73"}
        titulo_panel = lambda valor: rf"$R_c={valor:g}$ mm"
        etiqueta_linea = lambda valor: rf"$\eta={valor:g}$"
        filas_heat = [0.0045, 0.0050, 0.0055]
        columnas_heat = [0.55, 0.67, 0.80]
        fila_key, col_key = "Rc_mm", "eta"
        etiqueta_y_heat = r"$R_c$ [mm]"
        etiqueta_x_heat = r"$\eta=R_n/R_c$"

    fig, ejes = plt.subplots(1, 3, figsize=(10.0, 3.25), sharex=True, sharey=True)
    for eje, panel in zip(ejes, paneles):
        casos_panel = sorted(
            {
                fila["caso"]
                for fila in trayectorias
                if math.isclose(_numero(fila[parametro_panel]), panel)
            }
        )
        for codigo in casos_panel:
            filas = [fila for fila in trayectorias if fila["caso"] == codigo]
            if not filas:
                continue
            valor_linea = _numero(filas[0][parametro_linea])
            color = colores[valor_linea]
            alpha = 1.0 if filas[0]["estado"] == "normal" else 0.55
            t = [_numero(fila["tiempo_ms"]) for fila in filas]
            lc = [_numero(fila["lambda_citoplasma"]) for fila in filas]
            ln = [_numero(fila["lambda_nucleo"]) for fila in filas]
            eje.plot(t, lc, color=color, lw=1.6, alpha=alpha)
            eje.plot(t, ln, color=color, lw=1.4, ls="--", alpha=alpha)
        eje.set_title(titulo_panel(panel))
        eje.set_xlabel("Tiempo [ms]")
        eje.grid(True, color="#d9d9d9", lw=0.55)
    ejes[0].set_ylabel(r"Aspiración normalizada, $\lambda=(x-x_0)/R_c$")

    leyenda_parametros = [
        Line2D([0], [0], color=colores[valor], lw=1.7, label=etiqueta_linea(valor))
        for valor in colores
    ]
    leyenda_componentes = [
        Line2D([0], [0], color="#333333", lw=1.7, label="Citoplasma"),
        Line2D([0], [0], color="#333333", lw=1.5, ls="--", label="Núcleo"),
    ]
    fig.legend(
        handles=leyenda_parametros + leyenda_componentes,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.04),
        ncol=len(leyenda_parametros) + len(leyenda_componentes),
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    base = destino / f"{campana}_trayectorias_normalizadas"
    fig.savefig(base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    matriz = np.full((len(filas_heat), len(columnas_heat)), np.nan)
    estados = [["pendiente" for _ in columnas_heat] for _ in filas_heat]
    codigos = [["" for _ in columnas_heat] for _ in filas_heat]
    for fila in resumen:
        fv = _numero(fila[fila_key])
        cv = _numero(fila[col_key])
        if math.isnan(fv) or math.isnan(cv):
            continue
        i = next(i for i, valor in enumerate(filas_heat) if math.isclose(valor, fv))
        j = next(
            j for j, valor in enumerate(columnas_heat) if math.isclose(valor, cv)
        )
        matriz[i, j] = _numero(fila["lambda_citoplasma_ref"])
        estados[i][j] = str(fila["estado"])
        codigos[i][j] = str(fila["caso"])

    fig, eje = plt.subplots(figsize=(5.2, 3.8))
    cmap = plt.colormaps["viridis"].copy()
    cmap.set_bad("#e6e6e6")
    imagen = eje.imshow(matriz, cmap=cmap, aspect="auto", origin="upper")
    eje.set_xticks(range(len(columnas_heat)), [f"{x:g}" for x in columnas_heat])
    eje.set_yticks(range(len(filas_heat)), [f"{x:g}" for x in filas_heat])
    eje.set_xlabel(etiqueta_x_heat)
    eje.set_ylabel(etiqueta_y_heat)
    tiempo_ref = next(
        (_numero(f["tiempo_referencia_ms"]) for f in resumen if f["estado"] != "pendiente"),
        math.nan,
    )
    eje.set_title(
        rf"Aspiración del citoplasma a $t={tiempo_ref:.3f}$ ms, "
        rf"$\lambda_c=(x_c-x_{{c,0}})/R_c$"
    )
    valores_validos = matriz[~np.isnan(matriz)]
    punto_medio = (
        (float(np.min(valores_validos)) + float(np.max(valores_validos))) / 2
        if valores_validos.size
        else 0
    )
    for i in range(len(filas_heat)):
        for j in range(len(columnas_heat)):
            if math.isnan(float(matriz[i, j])):
                texto, color = "pendiente", "#555555"
            else:
                marca = "" if estados[i][j] == "normal" else "†"
                texto = f"{matriz[i, j]:.2f}{marca}"
                color = "white" if matriz[i, j] < punto_medio else "black"
            eje.text(j, i, texto, ha="center", va="center", color=color, fontsize=8)
            eje.text(
                j,
                i + 0.32,
                codigos[i][j],
                ha="center",
                va="center",
                color=color,
                fontsize=6.5,
            )
    barra = fig.colorbar(imagen, ax=eje)
    barra.set_label(r"$\lambda_c$")
    if any(
        fila["estado"] not in ("normal", "pendiente")
        for fila in resumen
    ):
        fig.text(
            0.01,
            0.01,
            "† Resultado utilizable hasta el instante comparado, pero sin terminación normal.",
            fontsize=7,
        )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    base = destino / f"{campana}_mapa_calor"
    fig.savefig(base.with_suffix(".png"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def argumentos() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--campana",
        choices=("material", "geometria", "todas"),
        default="todas",
        help="Campaña que se procesará.",
    )
    parser.add_argument(
        "--tiempo-referencia-ms",
        type=float,
        default=0.035,
        help="Instante común usado en los mapas de calor (por defecto: 0.035 ms).",
    )
    parser.add_argument(
        "--lspp",
        type=Path,
        default=LSPP_PREDETERMINADO,
        help="Ruta al ejecutable de LS-PrePost.",
    )
    parser.add_argument(
        "--solo-tablas",
        action="store_true",
        help="No genera PNG/PDF; útil si matplotlib no está instalado.",
    )
    return parser.parse_args()


def main() -> int:
    args = argumentos()
    if not args.lspp.exists():
        print(f"No existe LS-PrePost en {args.lspp}", file=sys.stderr)
        return 2
    SALIDA.mkdir(parents=True, exist_ok=True)
    seleccion = (
        [("material", CASOS_MATERIAL), ("geometria", CASOS_GEOMETRIA)]
        if args.campana == "todas"
        else [
            (
                args.campana,
                CASOS_MATERIAL if args.campana == "material" else CASOS_GEOMETRIA,
            )
        ]
    )
    tiempo_ref_s = args.tiempo_referencia_ms / 1000
    for campana, casos in seleccion:
        print(f"\nProcesando campaña: {campana}", flush=True)
        resumen, trayectorias = procesar_campana(
            campana, casos, args.lspp, tiempo_ref_s
        )
        if not args.solo_tablas:
            generar_figuras(campana, resumen, trayectorias)
    print(f"\nResultados escritos en: {SALIDA}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
