#!/usr/bin/env python3
"""Grafica el cambio de separación axial entre núcleo y citoplasma.

La entrada es ``trayectorias_material.csv``, generada por
``analizar_campanas.py``. Para cada caso se representa

    Delta d_nc / Rc = [(x_c - x_n) - (x_c,0 - x_n,0)] / Rc

frente al avance normalizado del citoplasma durante su tránsito por la
constricción. No se aplica suavizado.
"""

from __future__ import annotations

import csv
import math

import matplotlib.pyplot as plt

from avance_materialidad import (
    COLORES,
    ESTILOS,
    SALIDA,
    X_ENTRADA_MM,
    X_SALIDA_MM,
    formato_decimal,
    leer_trayectorias,
    numero,
)


def interpolar_cruce(
    filas: list[dict[str, str]], posicion_objetivo: float
) -> tuple[float, float]:
    """Interpola el tiempo y Delta d_nc/Rc en un cruce del frente."""
    for fila_0, fila_1 in zip(filas, filas[1:]):
        x_0 = numero(fila_0, "x_citoplasma_mm")
        x_1 = numero(fila_1, "x_citoplasma_mm")
        if x_0 <= posicion_objetivo <= x_1 and not math.isclose(x_0, x_1):
            fraccion = (posicion_objetivo - x_0) / (x_1 - x_0)
            t_0 = numero(fila_0, "tiempo_ms")
            t_1 = numero(fila_1, "tiempo_ms")
            s_0 = numero(fila_0, "rezago_nucleo")
            s_1 = numero(fila_1, "rezago_nucleo")
            tiempo = t_0 + fraccion * (t_1 - t_0)
            separacion = s_0 + fraccion * (s_1 - s_0)
            return tiempo, separacion
    raise ValueError(
        f"La trayectoria no cruza x={posicion_objetivo:.4f} mm."
    )


def construir_curva(
    filas: list[dict[str, str]],
) -> tuple[list[float], list[float]]:
    """Construye xi_c y Delta d_nc/Rc dentro de la constricción."""
    filas = sorted(filas, key=lambda fila: numero(fila, "tiempo_ms"))
    t_entrada, s_entrada = interpolar_cruce(filas, X_ENTRADA_MM)
    t_salida, s_salida = interpolar_cruce(filas, X_SALIDA_MM)
    radio = numero(filas[0], "Rc_mm")

    avance = [0.0]
    separacion = [s_entrada]
    for fila in filas:
        tiempo = numero(fila, "tiempo_ms")
        if t_entrada < tiempo < t_salida:
            avance.append(
                (numero(fila, "x_citoplasma_mm") - X_ENTRADA_MM) / radio
            )
            separacion.append(numero(fila, "rezago_nucleo"))

    avance.append((X_SALIDA_MM - X_ENTRADA_MM) / radio)
    separacion.append(s_salida)
    return avance, separacion


def generar_figuras() -> None:
    trayectorias = leer_trayectorias()
    SALIDA.mkdir(parents=True, exist_ok=True)

    curvas: dict[str, tuple[list[float], list[float]]] = {}
    for codigo, filas in trayectorias.items():
        curvas[codigo] = construir_curva(filas)

    valores = [valor for _, separacion in curvas.values() for valor in separacion]
    amplitud = max(valores) - min(valores)
    margen = 0.08 * amplitud if amplitud else 0.01
    limite_inferior = min(0.0, min(valores) - margen)
    limite_superior = max(0.0, max(valores) + margen)

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
    for en_kpa in (100.0, 225.0, 300.0):
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
            avance, separacion = curvas[codigo]
            etiqueta = (
                rf"$E_c={formato_decimal(ec_kpa)}\ \mathrm{{kPa}},"
                rf"\ \Gamma_E={formato_decimal(gamma)}$"
            )
            eje.plot(
                avance,
                separacion,
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
                    "separacion_entrada": separacion[0],
                    "separacion_minima": min(separacion),
                    "separacion_maxima": max(separacion),
                    "separacion_salida": separacion[-1],
                    "estado": filas[0]["estado"],
                }
            )

        eje.axhline(0.0, color="#777777", linewidth=0.8, linestyle="--")
        eje.set_xlabel(r"$\xi_c$ [-]")
        eje.set_ylabel(r"$\Delta d_{nc}/R_c$ [-]")
        eje.set_title(rf"$E_n={en_kpa:g}\ \mathrm{{kPa}}$")
        eje.set_xlim(0.0, (X_SALIDA_MM - X_ENTRADA_MM) / 0.005)
        eje.set_ylim(limite_inferior, limite_superior)
        eje.grid(True, color="#d9d9d9", linewidth=0.5)
        eje.spines["top"].set_visible(False)
        eje.spines["right"].set_visible(False)
        eje.legend(frameon=False, loc="best")
        fig.tight_layout()

        base = SALIDA / f"separacion_interna_En{en_kpa:g}"
        fig.savefig(base.with_suffix(".png"), bbox_inches="tight")
        fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)

    campos = [
        "caso",
        "En_kPa",
        "Ec_kPa",
        "Gamma",
        "separacion_entrada",
        "separacion_minima",
        "separacion_maxima",
        "separacion_salida",
        "estado",
    ]
    with (SALIDA.parent / "resumen_separacion_interna.csv").open(
        "w", encoding="utf-8", newline=""
    ) as archivo:
        escritor = csv.DictWriter(archivo, fieldnames=campos)
        escritor.writeheader()
        escritor.writerows(resumen)


if __name__ == "__main__":
    generar_figuras()
