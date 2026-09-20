"""Resultados VE-01--VE-09: extrae, comprueba y grafica los archivos existentes.

Uso desde cualquier directorio: py resultados_viscoelasticidad.py
Dependencias: numpy, matplotlib y LS-PrePost (solo para extraer d3plot).
No ejecuta LS-DYNA ni modifica las simulaciones. --extraer renueva la caché.
Las series originales se conservan. Los cruces se interpolan linealmente y los
casos sin salida no se extrapolan. Las curvas mantienen su orden temporal.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

import analizar_campanas as ac

LS = Path(__file__).resolve().parents[2]
OUT = LS / "resultados_procesados" / "viscoelasticidad"
LATEX = LS.parent / "tesis_latex"
RC, XENT, XSAL = 0.005, 0.0265, 0.046
# Referencia NC-08, obtenida con los mismos planos y la misma interpolación.
TC_REF_MS = 0.028601481093
COLORES = ("#0072B2", "#E69F00", "#009E73")
ESTILOS = ("-", "--", "-.")


def materiales(modelo):
    bloques, actual = [], None
    with modelo.open(encoding="latin-1") as f:
        for linea in f:
            if linea.startswith("*"):
                actual = [] if linea.strip() == "*MAT_SOFT_TISSUE_VISCO_TITLE" else None
                if actual is not None:
                    bloques.append(actual)
            elif actual is not None and linea.strip() and not linea.lstrip().startswith("$"):
                actual.append(linea.rstrip())
    def campos(linea):
        return [float(linea[i:i+10].strip() or 0) for i in range(0, 60, 10)]
    resultado = {}
    for b in bloques:
        resultado[int(b[1][:10])] = {
            "C1_MPa": float(b[1][20:30]), "XK_MPa": float(b[2][:10]),
            "S": campos(b[5]), "T_s": campos(b[6]),
        }
    assert set(resultado) == {2, 5}, "Se requieren los dos materiales viscoelásticos"
    c, n = resultado[2], resultado[5]
    assert c["S"] == n["S"] and c["T_s"] == n["T_s"]
    assert math.isclose(c["C1_MPa"], .025) and math.isclose(n["C1_MPa"], .05)
    assert c["XK_MPa"] == 2500 and n["XK_MPa"] == 5000
    assert c["T_s"][1] == 1e12 and all(s == 0 for s in c["S"][2:])
    assert math.isclose(sum(c["S"]), 1)
    return resultado


def huella(modelo, carpeta):
    h = hashlib.sha256()
    with modelo.open("rb") as f:
        for trozo in iter(lambda: f.read(1024*1024), b""):
            h.update(trozo)
    return {"modelo_sha256": h.hexdigest(), "d3plot": [
        [p.name, p.stat().st_size, p.stat().st_mtime_ns]
        for p in sorted(carpeta.glob("d3plot*")) if p.is_file()]}


def extraer(caso, carpeta, modelo, nodos, renovar):
    destino = OUT / "series"
    destino.mkdir(parents=True, exist_ok=True)
    rutas = {c: destino / f"{caso}_{c}_x.csv" for c in nodos}
    registro = destino / f"{caso}_fuentes.json"
    firma = huella(modelo, carpeta)
    vigente = registro.exists() and json.loads(registro.read_text()) == firma
    if renovar or not vigente or not all(p.exists() for p in rutas.values()):
        with tempfile.TemporaryDirectory(prefix="extraer_visco_", dir=OUT) as temporal:
            trabajo = Path(temporal)
            for p in [modelo, *sorted(carpeta.glob("d3plot*"))]:
                os.link(p, trabajo / p.name)
            archivo = trabajo / "extraer.cfile"
            salidas = {c: trabajo / f"{c}.csv" for c in nodos}
            ac.escribir_cfile(archivo, modelo, nodos, salidas)
            inicio = subprocess.STARTUPINFO()
            inicio.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            inicio.wShowWindow = 0
            subprocess.run([str(ac.LSPP_PREDETERMINADO), f"c={archivo}"],
                           cwd=trabajo, startupinfo=inicio, timeout=300, check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for c, ruta in salidas.items():
                ac.leer_csv_lspp(ruta)  # validar antes de reemplazar la caché
                shutil.copy2(ruta, rutas[c])
        registro.write_text(json.dumps(firma, indent=2), encoding="utf-8")
    c, n = (np.asarray(ac.leer_csv_lspp(rutas[k])) for k in ("citoplasma", "nucleo"))
    assert c.shape == n.shape and np.array_equal(c[:, 0], n[:, 0])
    assert np.isfinite(c).all() and np.isfinite(n).all()
    assert (np.diff(c[:, 0]) > 0).all()
    assert abs(c[0, 1]-nodos["citoplasma"][1]) < 1e-8
    assert abs(n[0, 1]-nodos["nucleo"][1]) < 1e-8
    assert abs(c[-1, 0]*1000-.15) < 1e-6
    return c[:, 0]*1000, c[:, 1], n[:, 1]


def cruce(t, x, objetivo):
    indices = np.flatnonzero((x[:-1] <= objetivo) & (x[1:] >= objetivo) & (np.diff(x) > 0))
    if len(indices) == 0:
        return None
    i = indices[0]
    return float(t[i]+(objetivo-x[i])/(x[i+1]-x[i])*(t[i+1]-t[i]))


def tramo(t, *variables, inicio, fin):
    # Inserta únicamente extremos interpolados, sin extrapolación ni filtros.
    assert t[0] <= inicio <= fin <= t[-1]
    tiempos = np.concatenate(([inicio], t[(t > inicio) & (t < fin)], [fin]))
    return (tiempos, *(np.interp(tiempos, t, v) for v in variables))


def energia(carpeta, inicio_ms, fin_ms):
    valores, tiempo = [], None
    with (carpeta / "d3hsp").open(encoding="latin-1") as f:
        for linea in f:
            if re.match(r"^\s*time\.{3,}", linea):
                tiempo = float(linea.split()[-1])*1000
            elif linea.strip().startswith("energy ratio w/o eroded energy."):
                if tiempo is not None and inicio_ms <= tiempo <= fin_ms:
                    valores.append(float(linea.split()[-1]))
    if not valores:
        raise ValueError(f"Sin registros energéticos en {carpeta}")
    return min(valores), max(valores), len(valores)


def decimal(x, cifras=5):
    return f"{x:.{cifras}f}".replace(".", "{,}")


def tabla(filas):
    texto = [r"% Generado por resultados_viscoelasticidad.py a partir de las series nodales.",
             r"\begin{tabular}{lcccc}", r"\toprule",
             r"Caso & $t_{\mathrm{ent}}$ [ms] & $\Delta t_{\mathrm{tr}}$ [ms] & Salida observada & $100\chi_f$ [\%] \\",
             r"\midrule"]
    for r in filas:
        duracion = r["duracion_ms"]
        if duracion is None:
            # Redondear hacia abajo conserva el significado de la cota inferior.
            cota = math.floor(r["cota_inferior_ms"]*1e5)/1e5
            d = "$>" + decimal(cota) + "$"
            progreso = "$"+decimal(r["progreso_final_pct"], 2)+"$"
        else:
            d, progreso = "$"+decimal(duracion)+"$", "---"
        texto.append(" & ".join([r["caso"], "$"+decimal(r["t_entrada_ms"])+"$", d,
                                  "Sí" if duracion is not None else "No", progreso])+r" \\")
    texto.extend([r"\bottomrule", r"\end{tabular}"])
    destino = LATEX / "tablas"
    destino.mkdir(exist_ok=True)
    (destino / "visco_resultados_filas.tex").write_text("\n".join(texto)+"\n", encoding="utf-8")
    minimo, maximo = min(r["RE_min"] for r in filas), max(r["RE_max"] for r in filas)
    fila = "Viscoelasticidad & "+decimal(minimo, 4).replace("{,}", ",")+" & "+decimal(maximo, 4).replace("{,}", ",")
    fila += " & "+decimal(max(abs(1-minimo), abs(1-maximo))*100, 2).replace("{,}", ",")+r" \\"+"\n"
    (destino / "visco_energia_fila.tex").write_text(
        r"\def\FilaViscoEnergia{"+fila.rstrip()+"}\n", encoding="utf-8")


def graficar(series):
    destino = OUT / "figuras"
    destino.mkdir(exist_ok=True)
    tesis_fig = LATEX / "figures" / "c4"
    tesis_fig.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "dejavuserif",
                         "font.size": 9, "axes.labelsize": 10, "legend.fontsize": 8,
                         "figure.dpi": 120, "savefig.dpi": 300, "pdf.fonttype": 42})
    for grupo in range(3):
        casos = series[grupo*3:grupo*3+3]
        s1, s2 = casos[0][0]["S1"], casos[0][0]["S2"]
        for tipo in ("avance", "separacion"):
            fig, ax = plt.subplots(figsize=(5.7, 3.5))
            for j, (r, t, xi, sep) in enumerate(casos):
                fin = r["t_salida_ms"] or float(t[-1])
                inicio = float(t[0]) if tipo == "avance" else r["t_entrada_ms"]
                tt, xx, ss = tramo(t, xi, sep, inicio=inicio, fin=fin)
                etiqueta = rf"$T_1={decimal(r['T1_ms'], 7).rstrip('0')}\ \mathrm{{ms}}$"
                ax.plot(tt if tipo == "avance" else xx, xx if tipo == "avance" else ss,
                        lw=1.8, color=COLORES[j], ls=ESTILOS[j], label=etiqueta)
            ax.set_title(rf"$S_1={decimal(s1, 2)},\quad S_2={decimal(s2, 2)}$", fontsize=10)
            if tipo == "avance":
                ax.set(xlabel=r"$t$ [ms]", ylabel=r"$\xi_c$ [-]", xlim=(0, .153), ylim=(-1.4, 4.05))
                ax.set_xticks(np.arange(0, .151, .03))
                ax.set_yticks([-1, 0, 1, 2, 3, 4])
                ax.legend(frameon=False, loc="upper right")
            else:
                ax.set(xlabel=r"$\xi_c$ [-]", ylabel=r"$\Delta d_{nc}/R_c$ [-]", xlim=(-.04, 3.95), ylim=(-.02, .032))
                ax.set_xticks([0, 1, 2, 3])
                ax.set_yticks([-.02, -.01, 0, .01, .02, .03])
                ax.legend(frameon=False, loc="lower left")
            for eje in (ax.xaxis, ax.yaxis):
                eje.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}".replace(".", ",")))
            ax.grid(color="#d9d9d9", lw=.5)
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            nombre = f"visco_{tipo}_S2_{s2:.2f}".replace(".", "p")
            for extension in ("png", "pdf"):
                ruta = destino / f"{nombre}.{extension}"
                fig.savefig(ruta, bbox_inches="tight")
                if extension == "pdf":
                    shutil.copy2(ruta, tesis_fig / ruta.name)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extraer", action="store_true", help="Reextraer las series con LS-PrePost")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    filas, series, trayectorias, fuentes = [], [], [], []
    for i in range(1, 10):
        caso, carpeta = f"VE-{i:02}", LS / "resultados" / "visco" / str(i)
        modelos = list(carpeta.glob("*.k"))
        assert len(modelos) == 1, f"Modelo ambiguo o ausente en {carpeta}"
        modelo = modelos[0]
        print(f"{caso}: lectura de resultados existentes", flush=True)
        mat = materiales(modelo)
        nodos = ac.nodos_frontales(modelo)
        t, xc, xn = extraer(caso, carpeta, modelo, nodos, args.extraer)
        tin, tout = cruce(t, xc, XENT), cruce(t, xc, XSAL)
        assert tin is not None, f"{caso}: ingreso no registrado"
        fin = tout if tout is not None else float(t[-1])
        sep = ((xc-xn)-(nodos["citoplasma"][1]-nodos["nucleo"][1]))/RC
        xi = (xc-XENT)/RC
        _, _, ss = tramo(t, xi, sep, inicio=tin, fin=fin)
        emin, emax, ne = energia(carpeta, tin, fin)
        m = mat[2]
        estado = ac.estado_simulacion(carpeta)
        assert estado == "normal", f"Revisar terminación de {caso}: {estado}"
        r = {"caso": caso, "S1": m["S"][0], "S2": m["S"][1],
             "T1_ms": m["T_s"][0]*1000, "De_ref": m["T_s"][0]*1000/TC_REF_MS,
             "t_entrada_ms": tin, "t_salida_ms": tout,
             "duracion_ms": tout-tin if tout is not None else None,
             "cota_inferior_ms": float(t[-1])-tin if tout is None else None,
             "t_final_ms": float(t[-1]), "salida_observada": tout is not None,
             "progreso_final_pct": float((xc[-1]-XENT)/(XSAL-XENT)*100) if tout is None else None,
             "separacion_min": float(min(ss)), "separacion_max": float(max(ss)),
             "separacion_fin_intervalo": float(ss[-1]), "RE_min": emin, "RE_max": emax,
             "registros_energia": ne, "estados_d3plot": len(t), "terminacion": estado}
        filas.append(r)
        series.append((r, t, xi, sep))
        fuentes.append({"caso": caso, "modelo": modelo.relative_to(LS).as_posix(),
                        "materiales": mat, "nodos": nodos})
        trayectorias.extend({"caso": caso, "t_ms": float(a), "xc_mm": float(b),
                            "xn_mm": float(c), "xi_c": float(d), "delta_d_sobre_Rc": float(e)}
                           for a,b,c,d,e in zip(t,xc,xn,xi,sep))
    assert len(filas) == 9
    for nombre, registros in (("resumen_viscoelasticidad", filas), ("trayectorias_viscoelasticidad", trayectorias)):
        with (OUT / f"{nombre}.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(registros[0]))
            w.writeheader()
            w.writerows(registros)
    (OUT / "resumen_viscoelasticidad.json").write_text(json.dumps(filas, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUT / "fuentes.json").write_text(json.dumps({"Rc_mm": RC, "x_ent_mm": XENT,
        "x_sal_mm": XSAL, "tc_ref_ms": TC_REF_MS, "casos": fuentes}, indent=2), encoding="utf-8")
    tabla(filas)
    graficar(series)
    print(f"Salida del frente observada en {sum(r['salida_observada'] for r in filas)} de 9 casos.")
    print(f"Seis figuras y tablas guardadas en {OUT}", flush=True)


if __name__ == "__main__":
    main()
