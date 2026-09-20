#!/usr/bin/env python3
"""Grafica el progreso del citoplasma dentro de la constricción.

La entrada es ``trayectorias_geometria.csv``, generada por
``analizar_campanas.py`` a partir de los archivos d3plot. Para cada caso se
representa

    chi = (x_c(t) - x_ent) / (x_sal - x_ent)

frente a ``Delta t_ent = t - t_ent``. Los cruces de entrada y salida se interpolan
linealmente entre estados consecutivos. No se aplica suavizado ni filtrado.
Si un caso no alcanza la salida, la curva termina en el último estado
disponible y no se extrapola.
"""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


TESIS_LSDYNA = Path(__file__).resolve().parents[2]
ENTRADA = (
    TESIS_LSDYNA
    / "resultados_procesados"
    / "geometria"
    / "trayectorias_geometria.csv"
)
DIRECTORIO_RESULTADOS = TESIS_LSDYNA / "resultados_procesados" / "geometria"
DIRECTORIO_FIGURAS = DIRECTORIO_RESULTADOS / "figuras"

X_ENTRADA_MM = 0.0265
X_SALIDA_MM = 0.0460
LONGITUD_CONSTRICCION_MM = X_SALIDA_MM - X_ENTRADA_MM

COLORES = {0.55: "#0072B2", 0.67: "#E69F00", 0.80: "#009E73"}
ESTILOS = {0.55: "-", 0.67: "--", 0.80: "-."}


def numero(fila: dict[str, str], campo: str) -> float:
    return float(fila[campo])


def primer_cruce(
    filas: list[dict[str, str]],
    posicion_objetivo: float,
    tiempo_minimo: float = -math.inf,
) -> tuple[float, int] | None:
    """Devuelve el primer cruce ascendente posterior a ``tiempo_minimo``."""
    for indice, (fila_0, fila_1) in enumerate(zip(filas, filas[1:])):
        t_0 = numero(fila_0, "tiempo_ms")
        t_1 = numero(fila_1, "tiempo_ms")
        if t_1 < tiempo_minimo:
            continue
        x_0 = numero(fila_0, "x_citoplasma_mm")
        x_1 = numero(fila_1, "x_citoplasma_mm")
        if x_0 <= posicion_objetivo <= x_1 and not math.isclose(x_0, x_1):
            fraccion = (posicion_objetivo - x_0) / (x_1 - x_0)
            tiempo = t_0 + fraccion * (t_1 - t_0)
            if tiempo >= tiempo_minimo:
                return tiempo, indice
    return None


def construir_curva(
    filas: list[dict[str, str]],
) -> tuple[list[float], list[float], dict[str, float | bool | str]]:
    """Construye la curva desde la entrada hasta la salida o el último estado."""
    filas = sorted(filas, key=lambda fila: numero(fila, "tiempo_ms"))
    entrada = primer_cruce(filas, X_ENTRADA_MM)
    if entrada is None:
        raise ValueError("La trayectoria no alcanza la entrada de la constricción.")

    t_entrada, _ = entrada
    salida = primer_cruce(filas, X_SALIDA_MM, tiempo_minimo=t_entrada)
    completa = salida is not None
    t_final = salida[0] if salida else numero(filas[-1], "tiempo_ms")

    tau = [0.0]
    progreso = [0.0]
    for fila in filas:
        tiempo = numero(fila, "tiempo_ms")
        if t_entrada < tiempo < t_final:
            tau.append(tiempo - t_entrada)
            progreso.append(
                (numero(fila, "x_citoplasma_mm") - X_ENTRADA_MM)
                / LONGITUD_CONSTRICCION_MM
            )

    if completa:
        progreso_final = 1.0
    else:
        progreso_final = (
            numero(filas[-1], "x_citoplasma_mm") - X_ENTRADA_MM
        ) / LONGITUD_CONSTRICCION_MM

    tau.append(t_final - t_entrada)
    progreso.append(progreso_final)
    metadatos: dict[str, float | bool | str] = {
        "t_entrada_ms": t_entrada,
        "t_salida_ms": salida[0] if salida else "",
        "duracion_constriccion_ms": salida[0] - t_entrada if salida else "",
        "progreso_observado": progreso_final,
        "completa_paso": completa,
    }
    return tau, progreso, metadatos


def leer_trayectorias() -> dict[str, list[dict[str, str]]]:
    if not ENTRADA.exists():
        raise FileNotFoundError(
            f"No existe {ENTRADA}. Ejecute primero analizar_campanas.py."
        )
    por_caso: dict[str, list[dict[str, str]]] = defaultdict(list)
    with ENTRADA.open("r", encoding="utf-8-sig", newline="") as archivo:
        for fila in csv.DictReader(archivo):
            por_caso[fila["caso"]].append(fila)
    return dict(por_caso)


def decimal_latex(valor: float, decimales: int) -> str:
    return f"{valor:.{decimales}f}".replace(".", "{,}")


def generar_figuras() -> None:
    trayectorias = leer_trayectorias()
    DIRECTORIO_FIGURAS.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "mathtext.fontset": "dejavuserif",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
        }
    )

    resumen: list[dict[str, str | float | bool]] = []
    for radio in (0.0045, 0.0050, 0.0055):
        fig, eje = plt.subplots(figsize=(5.7, 4.0))
        casos = [
            (codigo, filas)
            for codigo, filas in trayectorias.items()
            if math.isclose(numero(filas[0], "Rc_mm"), radio)
        ]
        casos.sort(key=lambda elemento: numero(elemento[1][0], "eta"))

        for codigo, filas in casos:
            eta = numero(filas[0], "eta")
            tau, progreso, datos = construir_curva(filas)
            eje.plot(
                tau,
                progreso,
                color=COLORES[eta],
                linestyle=ESTILOS[eta],
                linewidth=1.8,
                label=rf"$\eta={decimal_latex(eta, 2)}$",
            )
            resumen.append(
                {
                    "caso": codigo,
                    "Rc_mm": radio,
                    "eta": eta,
                    **datos,
                    "estado": filas[0]["estado"],
                }
            )

        eje.axhline(1.0, color="#777777", linewidth=0.8, linestyle="--")
        eje.set_xlabel(r"$\Delta t_{\mathrm{ent}}$ [ms]")
        eje.set_ylabel(r"$\chi$ [-]")
        eje.set_title(
            rf"$R_c={decimal_latex(radio, 4)}\ \mathrm{{mm}}$"
        )
        eje.set_xlim(0.0, 0.050)
        eje.set_ylim(0.0, 1.04)
        eje.grid(True, color="#d9d9d9", linewidth=0.5)
        eje.spines["top"].set_visible(False)
        eje.spines["right"].set_visible(False)
        eje.legend(frameon=False, loc="best")
        fig.tight_layout()

        sufijo = f"{radio:.4f}".replace(".", "p")
        base = DIRECTORIO_FIGURAS / f"progreso_constriccion_Rc{sufijo}"
        fig.savefig(base.with_suffix(".png"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)

    campos = [
        "caso",
        "Rc_mm",
        "eta",
        "t_entrada_ms",
        "t_salida_ms",
        "duracion_constriccion_ms",
        "progreso_observado",
        "completa_paso",
        "estado",
    ]
    with (DIRECTORIO_RESULTADOS / "resumen_progreso_constriccion.csv").open(
        "w", encoding="utf-8", newline=""
    ) as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=campos)
        escritor.writeheader()
        escritor.writerows(resumen)


if __name__ == "__main__":
    generar_figuras()
