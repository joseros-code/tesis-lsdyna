"""Genera la campaña viscoelástica VE-01 a VE-09.

La campaña reutiliza el caso hiperelástico NC-08 (En=300 kPa,
Ec=150 kPa) y reemplaza las tarjetas de núcleo y citoplasma por
*MAT_SOFT_TISSUE_VISCO. El tiempo característico se obtiene del resumen
de avance de NC-08 y los tiempos de relajación se definen mediante
De = T1/tc.

Ejemplos (desde tesis_lsdyna):
  py scripts/campaigns/barrido_visco.py --dry-run
  py scripts/campaigns/barrido_visco.py --only VE-05
  py scripts/campaigns/barrido_visco.py
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path


SCRIPTS_ROOT = Path(__file__).resolve().parent.parent
BASE = SCRIPTS_ROOT.parent

MODELO_BASE_PREDETERMINADO = (
    BASE
    / "generados"
    / "barrido_contraste_D18"
    / "modelo_NC-08_En300_Ec150_Gamma2p00.k"
)
RESUMEN_TC_PREDETERMINADO = (
    BASE
    / "resultados_procesados"
    / "material"
    / "resumen_avance_constriccion.csv"
)
SALIDA_PREDETERMINADA = BASE / "generados" / "barrido_viscoelastico_D18"

CASO_REFERENCIA = "NC-08"
EN_KPA = 300.0
EC_KPA = 150.0
GAMMA = EN_KPA / EC_KPA
ENDTIM_MS_PREDETERMINADO = 0.15
DT2MS_PREDETERMINADO = -1.11e-10
T2_S_PREDETERMINADO = 1.0e12

NIVELES_GINF = (0.50, 0.65, 0.80)
NIVELES_DE = (0.1, 1.0, 10.0)


def leer_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Genera nueve modelos viscoelásticos a partir de NC-08 "
            "mediante un barrido de g_inf y De."
        )
    )
    parser.add_argument(
        "--only",
        nargs="+",
        default=[],
        help="Genera solo los casos indicados, por ejemplo: --only VE-05 VE-09.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Valida y muestra la matriz sin escribir modelos .k.",
    )
    parser.add_argument(
        "--modelo-base",
        default=str(MODELO_BASE_PREDETERMINADO),
        help="Archivo hiperelástico NC-08 utilizado como plantilla.",
    )
    parser.add_argument(
        "--resumen-tc",
        default=str(RESUMEN_TC_PREDETERMINADO),
        help="CSV de avance que contiene los tiempos de NC-08.",
    )
    parser.add_argument(
        "--tc-ms",
        type=float,
        default=None,
        help="Sobrescribe tc en ms. Si se omite, se obtiene del CSV de NC-08.",
    )
    parser.add_argument(
        "--out-dir",
        default=str(SALIDA_PREDETERMINADA),
        help="Carpeta de salida de los modelos y archivos de trazabilidad.",
    )
    parser.add_argument(
        "--endtim-ms",
        type=float,
        default=ENDTIM_MS_PREDETERMINADO,
        help="Tiempo final de cada simulación en ms.",
    )
    parser.add_argument(
        "--dt2ms",
        type=float,
        default=DT2MS_PREDETERMINADO,
        help="DT2MS de *CONTROL_TIMESTEP. Use 0 para desactivar mass scaling.",
    )
    parser.add_argument(
        "--t2-s",
        type=float,
        default=T2_S_PREDETERMINADO,
        help="Tiempo del término casi permanente de la serie de Prony, en s.",
    )
    return parser.parse_args()


def resolver_ruta(valor: str) -> Path:
    ruta = Path(valor)
    if not ruta.is_absolute():
        ruta = BASE / ruta
    return ruta.resolve()


def obtener_tc(ruta_csv: Path, tc_ms_manual: float | None) -> dict[str, float | str]:
    if tc_ms_manual is not None:
        if tc_ms_manual <= 0.0:
            raise ValueError("--tc-ms debe ser positivo")
        return {
            "fuente": "valor_manual",
            "t_entrada_ms": math.nan,
            "t_salida_ms": math.nan,
            "tc_ms": tc_ms_manual,
        }

    if not ruta_csv.exists():
        raise FileNotFoundError(f"No existe el resumen para obtener tc: {ruta_csv}")

    with ruta_csv.open(encoding="utf-8", newline="") as archivo:
        for fila in csv.DictReader(archivo):
            if fila.get("caso", "").strip().upper() != CASO_REFERENCIA:
                continue
            en = float(fila["En_kPa"])
            ec = float(fila["Ec_kPa"])
            if not math.isclose(en, EN_KPA, abs_tol=1e-6):
                raise ValueError(f"{CASO_REFERENCIA}: En inesperado en el resumen: {en}")
            if not math.isclose(ec, EC_KPA, abs_tol=1e-6):
                raise ValueError(f"{CASO_REFERENCIA}: Ec inesperado en el resumen: {ec}")
            t_entrada = float(fila["t_entrada_ms"])
            t_salida = float(fila["t_salida_ms"])
            tc_reportado = float(fila["duracion_constriccion_ms"])
            tc_calculado = t_salida - t_entrada
            if not math.isclose(
                tc_reportado, tc_calculado, rel_tol=0.0, abs_tol=1e-12
            ):
                raise ValueError(
                    "El tc reportado no coincide con t_salida - t_entrada: "
                    f"{tc_reportado} != {tc_calculado}"
                )
            if tc_calculado <= 0.0:
                raise ValueError("El tiempo característico debe ser positivo")
            return {
                "fuente": str(ruta_csv),
                "t_entrada_ms": t_entrada,
                "t_salida_ms": t_salida,
                "tc_ms": tc_calculado,
            }

    raise ValueError(f"No se encontró {CASO_REFERENCIA} en {ruta_csv}")


def construir_casos(tc_ms: float) -> list[dict[str, float | int | str]]:
    casos: list[dict[str, float | int | str]] = []
    numero = 1
    for g_inf in NIVELES_GINF:
        for de in NIVELES_DE:
            t1_s = de * tc_ms / 1000.0
            casos.append(
                {
                    "case_id": f"VE-{numero:02d}",
                    "g_inf": g_inf,
                    "S1": 1.0 - g_inf,
                    "S2": g_inf,
                    "De": de,
                    "T1_s": t1_s,
                }
            )
            numero += 1
    return casos


def seleccionar_casos(
    casos: list[dict[str, float | int | str]], pedidos: list[str]
) -> list[dict[str, float | int | str]]:
    if not pedidos:
        return casos
    solicitados = {valor.upper() for valor in pedidos}
    conocidos = {str(caso["case_id"]) for caso in casos}
    desconocidos = sorted(solicitados - conocidos)
    if desconocidos:
        raise ValueError("Casos no reconocidos: " + ", ".join(desconocidos))
    return [caso for caso in casos if caso["case_id"] in solicitados]


def campo_float(valor: float, decimales: int = 4) -> str:
    texto = f"{valor:10.{decimales}E}"
    if len(texto) != 10:
        raise ValueError(f"El valor {valor} no cabe en un campo de 10 caracteres")
    return texto


def linea_seis_campos(valores: tuple[float, ...], cientifico: bool) -> str:
    if len(valores) != 6:
        raise ValueError("Se requieren exactamente seis valores")
    if cientifico:
        campos = [campo_float(valor) for valor in valores]
    else:
        campos = [f"{valor:10.6f}" for valor in valores]
    return "".join(campos) + "\n"


def siguiente_keyword(lineas: list[str], inicio: int) -> int:
    for indice in range(inicio, len(lineas)):
        if lineas[indice].lstrip().startswith("*"):
            return indice
    return len(lineas)


def datos_material(lineas: list[str], inicio: int, fin: int) -> list[int]:
    titulo_consumido = False
    indices: list[int] = []
    for indice in range(inicio + 1, fin):
        texto = lineas[indice].strip()
        if not texto or texto.startswith("$"):
            continue
        if not titulo_consumido:
            titulo_consumido = True
            continue
        indices.append(indice)
    return indices


def convertir_materiales(
    ruta: Path, caso: dict[str, float | int | str], t2_s: float
) -> None:
    lineas = ruta.read_text(encoding="utf-8", errors="ignore").splitlines(True)
    indices = [
        i
        for i, linea in enumerate(lineas)
        if linea.strip().upper() == "*MAT_SOFT_TISSUE_TITLE"
    ]
    if len(indices) != 2:
        raise ValueError(
            f"{ruta.name}: se esperaban 2 tarjetas *MAT_SOFT_TISSUE_TITLE; "
            f"se encontraron {len(indices)}"
        )

    mids: set[int] = set()
    for inicio in reversed(indices):
        fin = siguiente_keyword(lineas, inicio + 1)
        datos = datos_material(lineas, inicio, fin)
        if len(datos) != 4:
            raise ValueError(
                f"{ruta.name}: la tarjeta material iniciada en la línea "
                f"{inicio + 1} contiene {len(datos)} filas de datos; se esperaban 4"
            )
        try:
            mid = int(lineas[datos[0]][:10])
        except ValueError as exc:
            raise ValueError(f"{ruta.name}: no se pudo leer el MID") from exc
        mids.add(mid)
        lineas[inicio] = "*MAT_SOFT_TISSUE_VISCO_TITLE\n"
        insercion = datos[-1] + 1
        bloque_prony = [
            "$#      s1        s2        s3        s4        s5        s6\n",
            linea_seis_campos(
                (
                    float(caso["S1"]),
                    float(caso["S2"]),
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ),
                cientifico=False,
            ),
            "$#      t1        t2        t3        t4        t5        t6\n",
            linea_seis_campos(
                (float(caso["T1_s"]), t2_s, 0.0, 0.0, 0.0, 0.0),
                cientifico=True,
            ),
        ]
        lineas[insercion:insercion] = bloque_prony

    if mids != {2, 5}:
        raise ValueError(f"{ruta.name}: se esperaban MID 2 y 5; se encontraron {mids}")
    ruta.write_text("".join(lineas), encoding="utf-8")


def fijar_primer_campo(ruta: Path, keyword: str, valor: float) -> None:
    lineas = ruta.read_text(encoding="utf-8", errors="ignore").splitlines(True)
    for inicio, linea in enumerate(lineas):
        if not linea.strip().upper().startswith(keyword.upper()):
            continue
        for indice in range(inicio + 1, len(lineas)):
            texto = lineas[indice].strip()
            if texto.startswith("*"):
                break
            if texto and not texto.startswith("$"):
                original = lineas[indice].rstrip("\r\n")
                lineas[indice] = campo_float(valor) + original[10:] + "\n"
                ruta.write_text("".join(lineas), encoding="utf-8")
                return
    raise ValueError(f"{ruta.name}: no se encontró una fila de datos para {keyword}")


def fijar_dt2ms(ruta: Path, valor: float) -> None:
    lineas = ruta.read_text(encoding="utf-8", errors="ignore").splitlines(True)
    for inicio, linea in enumerate(lineas):
        if not linea.strip().upper().startswith("*CONTROL_TIMESTEP"):
            continue
        filas_datos: list[int] = []
        for indice in range(inicio + 1, len(lineas)):
            texto = lineas[indice].strip()
            if texto.startswith("*"):
                break
            if texto and not texto.startswith("$"):
                filas_datos.append(indice)
        if not filas_datos:
            break

        primera = lineas[filas_datos[0]].rstrip("\r\n").ljust(80)
        campos = [primera[i : i + 10] for i in range(0, 80, 10)]
        campos[4] = f"{valor:10.3E}"
        lineas[filas_datos[0]] = "".join(campos) + "\n"

        ruta.write_text("".join(lineas), encoding="utf-8")
        return
    raise ValueError(f"{ruta.name}: no se encontró *CONTROL_TIMESTEP")


def leer_fila_keyword(ruta: Path, keyword: str) -> str:
    lineas = ruta.read_text(encoding="utf-8", errors="ignore").splitlines()
    for inicio, linea in enumerate(lineas):
        if not linea.strip().upper().startswith(keyword.upper()):
            continue
        for fila in lineas[inicio + 1 :]:
            texto = fila.strip()
            if texto.startswith("*"):
                break
            if texto and not texto.startswith("$"):
                return fila
    raise ValueError(f"{ruta.name}: no se encontró una fila para {keyword}")


def valores_fijos(linea: str, cantidad: int = 6) -> list[float]:
    valores = []
    for indice in range(cantidad):
        campo = linea[indice * 10 : (indice + 1) * 10].strip()
        valores.append(float(campo) if campo else 0.0)
    return valores


def validar_modelo(
    ruta: Path,
    caso: dict[str, float | int | str],
    endtim_s: float,
    dt2ms: float,
    t2_s: float,
) -> None:
    texto = ruta.read_text(encoding="utf-8", errors="ignore")
    if texto.upper().count("*MAT_SOFT_TISSUE_VISCO_TITLE") != 2:
        raise ValueError(f"{ruta.name}: número incorrecto de materiales viscoelásticos")
    if "*MAT_SOFT_TISSUE_TITLE" in texto.upper().replace(
        "*MAT_SOFT_TISSUE_VISCO_TITLE", ""
    ):
        raise ValueError(f"{ruta.name}: quedó una tarjeta hiperelástica sin convertir")

    endtim = float(leer_fila_keyword(ruta, "*CONTROL_TERMINATION")[:10])
    timestep = valores_fijos(leer_fila_keyword(ruta, "*CONTROL_TIMESTEP"), 5)
    if not math.isclose(endtim, endtim_s, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError(f"{ruta.name}: ENDTIM={endtim}; esperado={endtim_s}")
    if not math.isclose(timestep[4], dt2ms, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError(f"{ruta.name}: DT2MS={timestep[4]}; esperado={dt2ms}")

    lineas = texto.splitlines()
    encontrados = 0
    for indice, linea in enumerate(lineas):
        if linea.strip().upper() != "*MAT_SOFT_TISSUE_VISCO_TITLE":
            continue
        fin = siguiente_keyword([valor + "\n" for valor in lineas], indice + 1)
        datos = datos_material([valor + "\n" for valor in lineas], indice, fin)
        if len(datos) != 6:
            raise ValueError(f"{ruta.name}: tarjeta viscoelástica incompleta")
        s = valores_fijos(lineas[datos[4]])
        t = valores_fijos(lineas[datos[5]])
        if not math.isclose(s[0], float(caso["S1"]), abs_tol=1e-8):
            raise ValueError(f"{ruta.name}: S1 incorrecto")
        if not math.isclose(s[1], float(caso["S2"]), abs_tol=1e-8):
            raise ValueError(f"{ruta.name}: S2 incorrecto")
        if not math.isclose(s[0] + s[1], 1.0, abs_tol=1e-8):
            raise ValueError(f"{ruta.name}: S1+S2 debe ser 1")
        if not math.isclose(t[0], float(caso["T1_s"]), rel_tol=5e-5):
            raise ValueError(f"{ruta.name}: T1 incorrecto: {t[0]}")
        if not math.isclose(t[1], t2_s, rel_tol=1e-8):
            raise ValueError(f"{ruta.name}: T2 incorrecto: {t[1]}")
        encontrados += 1
    if encontrados != 2:
        raise ValueError(f"{ruta.name}: se validaron {encontrados} materiales")


def tag_decimal(valor: float) -> str:
    return f"{valor:g}".replace(".", "p")


def nombre_modelo(caso: dict[str, float | int | str]) -> str:
    g_tag = f"{float(caso['g_inf']):.2f}".replace(".", "p")
    de_tag = tag_decimal(float(caso["De"]))
    return (
        f"modelo_{caso['case_id']}_En300_Ec150_"
        f"ginf{g_tag}_De{de_tag}.k"
    )


def escribir_trazabilidad(
    out_dir: Path,
    casos: list[dict[str, float | int | str]],
    tc: dict[str, float | str],
    modelo_base: Path,
    endtim_ms: float,
    dt2ms: float,
    t2_s: float,
) -> tuple[Path, Path, Path]:
    csv_path = out_dir / "casos_viscoelasticos.csv"
    campos = [
        "case_id",
        "En_kPa",
        "Ec_kPa",
        "Gamma_En_Ec",
        "g_inf",
        "S1",
        "S2",
        "De",
        "T1_s",
        "T1_ms",
        "T2_s",
        "tc_ms",
        "endtim_ms",
        "dt2ms",
        "archivo_modelo",
        "estado",
    ]
    with csv_path.open("w", encoding="utf-8", newline="") as archivo:
        writer = csv.DictWriter(archivo, fieldnames=campos)
        writer.writeheader()
        for caso in casos:
            writer.writerow(
                {
                    "case_id": caso["case_id"],
                    "En_kPa": f"{EN_KPA:.6f}",
                    "Ec_kPa": f"{EC_KPA:.6f}",
                    "Gamma_En_Ec": f"{GAMMA:.6f}",
                    "g_inf": f"{float(caso['g_inf']):.6f}",
                    "S1": f"{float(caso['S1']):.6f}",
                    "S2": f"{float(caso['S2']):.6f}",
                    "De": f"{float(caso['De']):.6f}",
                    "T1_s": f"{float(caso['T1_s']):.9E}",
                    "T1_ms": f"{float(caso['T1_s']) * 1000.0:.9E}",
                    "T2_s": f"{t2_s:.6E}",
                    "tc_ms": f"{float(tc['tc_ms']):.12f}",
                    "endtim_ms": f"{endtim_ms:.6f}",
                    "dt2ms": f"{dt2ms:.6E}",
                    "archivo_modelo": caso["archivo_modelo"],
                    "estado": caso["estado"],
                }
            )

    matriz_path = out_dir / "matriz_viscoelastica.txt"
    lineas = [
        "Campaña viscoelástica adimensional VE-01 a VE-09",
        "",
        f"Modelo hiperelástico base: {modelo_base}",
        f"Fuente de tc: {tc['fuente']}",
        f"t_entrada = {float(tc['t_entrada_ms']):.12f} ms",
        f"t_salida  = {float(tc['t_salida_ms']):.12f} ms",
        f"tc        = {float(tc['tc_ms']):.12f} ms",
        f"ENDTIM    = {endtim_ms:.6f} ms",
        f"DT2MS     = {dt2ms:.6E} s",
        f"T2        = {t2_s:.6E} s",
        "",
        "Los T1 son tiempos efectivos de la escala numérica y no tiempos",
        "de relajación biológicos medidos. La campaña evalúa sensibilidad",
        "respecto de De=T1/tc.",
        "",
        "Caso   g_inf   S1      S2      De       T1 [s]       Estado",
        "------  ------  ------  ------  -------  -----------  --------------------",
    ]
    for caso in casos:
        lineas.append(
            f"{str(caso['case_id']):<6}  {float(caso['g_inf']):>6.2f}  "
            f"{float(caso['S1']):>6.2f}  {float(caso['S2']):>6.2f}  "
            f"{float(caso['De']):>7.1f}  {float(caso['T1_s']):>11.5E}  "
            f"{caso['estado']}"
        )
    matriz_path.write_text("\n".join(lineas) + "\n", encoding="utf-8")

    orden_path = out_dir / "orden_simulacion.txt"
    orden = ["VE-05", "VE-04", "VE-06", "VE-02", "VE-08", "VE-01", "VE-03", "VE-07", "VE-09"]
    texto_orden = [
        "Orden recomendado de simulación",
        "",
        "Simular primero VE-05 (g_inf=0,65; De=1) como caso central.",
        "Revisar estabilidad, energía, deformación y masa añadida antes",
        "de continuar con los ocho casos restantes.",
        "",
    ]
    por_id = {str(caso["case_id"]): caso for caso in casos}
    for prioridad, case_id in enumerate(orden, start=1):
        caso = por_id[case_id]
        texto_orden.append(
            f"{prioridad}. {case_id}: g_inf={float(caso['g_inf']):.2f}, "
            f"De={float(caso['De']):.1f}, T1={float(caso['T1_s']):.5E} s"
        )
    orden_path.write_text("\n".join(texto_orden) + "\n", encoding="utf-8")
    return csv_path, matriz_path, orden_path


def main() -> None:
    args = leer_args()
    modelo_base = resolver_ruta(args.modelo_base)
    resumen_tc = resolver_ruta(args.resumen_tc)
    out_dir = resolver_ruta(args.out_dir)

    if not modelo_base.exists():
        raise FileNotFoundError(f"No existe el modelo base: {modelo_base}")
    if args.endtim_ms <= 0.0:
        raise ValueError("--endtim-ms debe ser positivo")
    if args.t2_s <= args.endtim_ms / 1000.0:
        raise ValueError("T2 debe ser mucho mayor que ENDTIM")

    tc = obtener_tc(resumen_tc, args.tc_ms)
    casos = construir_casos(float(tc["tc_ms"]))
    seleccionados = seleccionar_casos(casos, args.only)
    seleccion_ids = {str(caso["case_id"]) for caso in seleccionados}
    endtim_s = args.endtim_ms / 1000.0

    out_dir.mkdir(parents=True, exist_ok=True)
    print("Campaña viscoelástica VE-01 a VE-09")
    print("Modelo base:", modelo_base)
    print("Salida:", out_dir)
    print("tc:", f"{float(tc['tc_ms']):.12f}", "ms")
    print("ENDTIM:", args.endtim_ms, "ms")
    print("DT2MS:", args.dt2ms)
    print("Casos seleccionados:", len(seleccionados))

    for caso in casos:
        caso["archivo_modelo"] = nombre_modelo(caso)
        caso["estado"] = "omitido_por_filtro"

    for caso in seleccionados:
        ruta_salida = out_dir / str(caso["archivo_modelo"])
        print(
            f"{caso['case_id']}: g_inf={float(caso['g_inf']):.2f}, "
            f"De={float(caso['De']):.1f}, T1={float(caso['T1_s']):.6E} s"
        )
        if args.dry_run:
            caso["estado"] = "dry-run"
            continue
        shutil.copy2(modelo_base, ruta_salida)
        convertir_materiales(ruta_salida, caso, args.t2_s)
        fijar_primer_campo(ruta_salida, "*CONTROL_TERMINATION", endtim_s)
        fijar_dt2ms(ruta_salida, args.dt2ms)
        validar_modelo(
            ruta_salida,
            caso,
            endtim_s=endtim_s,
            dt2ms=args.dt2ms,
            t2_s=args.t2_s,
        )
        caso["estado"] = "generado_verificado"

    for caso in casos:
        if caso["case_id"] in seleccion_ids and caso["estado"] == "omitido_por_filtro":
            raise RuntimeError(f"{caso['case_id']}: no fue procesado")

    trazabilidad = escribir_trazabilidad(
        out_dir,
        casos,
        tc,
        modelo_base,
        endtim_ms=args.endtim_ms,
        dt2ms=args.dt2ms,
        t2_s=args.t2_s,
    )
    print("Archivos de trazabilidad:")
    for ruta in trazabilidad:
        print(" ", ruta)
    print("Campaña terminada.")


if __name__ == "__main__":
    main()
