#!/usr/bin/env python3
"""Grafica el avance normalizado durante el paso por la constricción.

La entrada es ``trayectorias_material.csv``, generada por
``analizar_campanas.py`` a partir de los archivos d3plot. Para cada caso se
interpola linealmente el instante en que el nodo frontal del citoplasma alcanza
la entrada y la salida de la constricción. No se aplica suavizado a la
trayectoria.
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
    / "material"
    / "trayectorias_material.csv"
)
SALIDA = TESIS_LSDYNA / "resultados_procesados" / "material" / "figuras"
X_ENTRADA_MM = 0.0265
X_SALIDA_MM = 0.0460

COLORES = {1.5: "#0072B2", 2.0: "#E69F00", 3.0: "#009E73"}
ESTILOS = {1.5: "-", 2.0: "--", 3.0: "-."}


def numero(fila: dict[str, str], campo: str) -> float:
    return float(fila[campo])


def interpolar_cruce(
    filas: list[dict[str, str]], posicion_objetivo: float
) -> tuple[float, int]:
    """Devuelve el primer tiempo de cruce y el índice del tramo que lo encierra."""
    for indice, (fila_0, fila_1) in enumerate(zip(filas, filas[1:])):
        x_0 = numero(fila_0, "x_citoplasma_mm")
        x_1 = numero(fila_1, "x_citoplasma_mm")
        if x_0 <= posicion_objetivo <= x_1 and not math.isclose(x_0, x_1):
            t_0 = numero(fila_0, "tiempo_ms")
            t_1 = numero(fila_1, "tiempo_ms")
            fraccion = (posicion_objetivo - x_0) / (x_1 - x_0)
            return t_0 + fraccion * (t_1 - t_0), indice
    raise ValueError(
        f"La trayectoria no cruza x={posicion_objetivo:.4f} mm."
    )


def construir_curva(
    filas: list[dict[str, str]],
) -> tuple[list[float], list[float], float, float]:
    """Construye Delta t_ent y xi_c entre ambos extremos de la constricción."""
    filas = sorted(filas, key=lambda fila: numero(fila, "tiempo_ms"))
    t_entrada, _ = interpolar_cruce(filas, X_ENTRADA_MM)
    t_salida, _ = interpolar_cruce(filas, X_SALIDA_MM)
    radio = numero(filas[0], "Rc_mm")

    tau = [0.0]
    avance = [0.0]
    for fila in filas:
        tiempo = numero(fila, "tiempo_ms")
        if t_entrada < tiempo < t_salida:
            tau.append(tiempo - t_entrada)
            avance.append(
                (numero(fila, "x_citoplasma_mm") - X_ENTRADA_MM) / radio
            )

    tau.append(t_salida - t_entrada)
    avance.append((X_SALIDA_MM - X_ENTRADA_MM) / radio)
    return tau, avance, t_entrada, t_salida


def formato_decimal(valor: float) -> str:
    if math.isclose(valor, round(valor)):
        return f"{valor:.0f}"
    if math.isclose(valor * 10, round(valor * 10)):
        texto = f"{valor:.1f}"
    else:
        texto = f"{valor:.2f}"
    return texto.replace(".", "{,}")


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


def generar_figuras() -> None:
    trayectorias = leer_trayectorias()
    SALIDA.mkdir(parents=True, exist_ok=True)

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

    resumen: list[dict[str, str | float]] = []
    grupos = (100.0, 225.0, 300.0)
    for en_kpa in grupos:
        fig, eje = plt.subplots(figsize=(5.7, 4.0))
        casos = [
            (codigo, filas)
            for codigo, filas in trayectorias.items()
            if math.isclose(numero(filas[0], "En_kPa"), en_kpa)
        ]
        casos.sort(key=lambda elemento: numero(elemento[1][0], "Gamma"))

        for codigo, filas in casos:
            gamma = numero(filas[0], "Gamma")
            ec_kpa = numero(filas[0], "Ec_kPa")
            tau, avance, t_entrada, t_salida = construir_curva(filas)
            etiqueta = (
                rf"$E_c={formato_decimal(ec_kpa)}\ \mathrm{{kPa}},"
                rf"\ \Gamma_E={formato_decimal(gamma)}$"
            )
            eje.plot(
                tau,
                avance,
                color=COLORES[gamma],
                linestyle=ESTILOS[gamma],
                linewidth=1.8,
                label=etiqueta,
            )
            resumen.append(
                {
                    "caso": codigo,
                    "En_kPa": en_kpa,
                    "Ec_kPa": ec_kpa,
                    "Gamma": gamma,
                    "t_entrada_ms": t_entrada,
                    "t_salida_ms": t_salida,
                    "duracion_constriccion_ms": t_salida - t_entrada,
                    "estado": filas[0]["estado"],
                }
            )

        eje.set_xlabel(r"$\Delta t_{\mathrm{ent}}$ [ms]")
        eje.set_ylabel(r"$\xi_c$ [-]")
        eje.set_title(rf"$E_n={en_kpa:g}\ \mathrm{{kPa}}$")
        eje.set_ylim(0.0, 4.05)
        eje.margins(x=0.015)
        eje.grid(True, color="#d9d9d9", linewidth=0.5)
        eje.spines["top"].set_visible(False)
        eje.spines["right"].set_visible(False)
        eje.legend(frameon=False, loc="best")
        fig.tight_layout()

        base = SALIDA / f"avance_normalizado_En{en_kpa:g}"
        fig.savefig(base.with_suffix(".png"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)

    campos = [
        "caso",
        "En_kPa",
        "Ec_kPa",
        "Gamma",
        "t_entrada_ms",
        "t_salida_ms",
        "duracion_constriccion_ms",
        "estado",
    ]
    with (SALIDA.parent / "resumen_avance_constriccion.csv").open(
        "w", encoding="utf-8", newline=""
    ) as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=campos)
        escritor.writeheader()
        escritor.writerows(resumen)


if __name__ == "__main__":
    generar_figuras()
