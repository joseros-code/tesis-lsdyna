# integrar_celula.py
# ------------------------------------------------------------
# Integra la celula nueva dentro del modelo base.
#
# Que modifica:
#   - elimina SOLO elementos solidos antiguos con PID=2
#   - agrega nodos nuevos de la celula generada
#   - agrega elementos nuevos PID=2 y PID=5
#   - agrega MID=5 y PART PID=5 si no existen
#
# Que mantiene intacto:
#   - fluido ALE PID=1
#   - shells/contacto PID=3 y PID=4
#   - nodos antiguos, sets, boxes, contacto, FSI, CONTROL, DATABASE
# ------------------------------------------------------------

import argparse
from pathlib import Path

from definir_materiales import (
    args_mat,
    leer_mat,
    mat_block,
)


# ============================================================
# 1. PARAMETROS
# ============================================================

base = Path(__file__).resolve().parents[2]


def resolver_ruta(ruta):
    ruta = Path(ruta)

    if ruta.is_absolute():
        return ruta

    return base / ruta


def leer_args():
    p = argparse.ArgumentParser(description="Agrega la celula nueva al modelo base.")
    p.add_argument(
        "--modelo",
        "--mod",
        dest="modelo",
        default="originales/MallaCelula_06.k",
        help="Archivo .k del modelo base.",
    )
    p.add_argument(
        "--celula",
        "--cell",
        dest="celula",
        default="generados/Celula_cuarto_esfera_Rc005_eta067.k",
        help="Archivo .k con la nueva celula.",
    )
    p.add_argument(
        "--salida",
        "--out",
        dest="salida",
        default="generados/MallaCelula_06_cuarto_esfera_nucleo.k",
        help="Archivo .k integrado de salida.",
    )
    p.add_argument(
        "--dt2ms",
        type=float,
        default=None,
        help="Valor de DT2MS para escalado de masa. Si se omite no se cambia.",
    )
    args_mat(p)

    return p.parse_args()


args = leer_args()

archivo_modelo_entrada = resolver_ruta(args.modelo)
archivo_celula_nueva = resolver_ruta(args.celula)
archivo_modelo_salida = resolver_ruta(args.salida)

pid_citoplasma = 2
pid_nucleo = 5

mid_nucleo = 5
mid_citoplasma = 2
secid_celula = 2

materiales = leer_mat(args)

ro_citoplasma = materiales["ro_citoplasma"]
c1_citoplasma = materiales["c1_citoplasma"]
xk_citoplasma = materiales["xk_citoplasma"]

ro_nucleo = materiales["ro_nucleo"]
c1_nucleo = materiales["c1_nucleo"]
xk_nucleo = materiales["xk_nucleo"]


# ============================================================
# 2. FUNCIONES DE TEXTO
# ============================================================

def separar_campos(linea):
    return linea.replace(",", " ").split()


def escribir_nodo(nid, coordenadas):
    x, y, z = coordenadas
    return f"{nid:8d}{x:16.8E}{y:16.8E}{z:16.8E}\n"


def escribir_elemento(eid, pid, conectividad):
    linea = f"{eid:8d}{pid:8d}"

    for nid in conectividad:
        linea += f"{nid:8d}"

    return linea + "\n"


def insertar_antes_de_keyword(lineas, bloque, keyword):
    posicion = len(lineas)

    for i, linea in enumerate(lineas):
        if linea.strip().upper().startswith(keyword.upper()):
            posicion = i
            break

    return lineas[:posicion] + bloque + lineas[posicion:]


def crear_bloque_part_nucleo():
    return [
        "$\n",
        "$ Parte agregada para el nucleo\n",
        "$\n",
        "*PART\n",
        "Nucleo_parametrico\n",
        "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n",
        f"{pid_nucleo:10d}{secid_celula:10d}{mid_nucleo:10d}"
        f"{0:10d}{0:10d}{0:10d}{0:10d}{0:10d}\n",
    ]


# ============================================================
# 3. LECTURA BASICA DE .K
# ============================================================

def leer_nodos_y_elementos_solidos(ruta):
    with open(ruta, "r", encoding="utf-8", errors="ignore") as f:
        lineas = f.readlines()

    nodos = {}
    elementos = {}
    max_nid = 0
    max_eid = 0

    leyendo_nodos = False
    leyendo_solidos = False

    for i, linea in enumerate(lineas):
        texto = linea.strip()

        if texto.upper().startswith("*NODE"):
            leyendo_nodos = True
            leyendo_solidos = False
            continue

        if texto.upper().startswith("*ELEMENT_SOLID"):
            leyendo_nodos = False
            leyendo_solidos = True
            continue

        if texto.startswith("*"):
            leyendo_nodos = False
            leyendo_solidos = False
            continue

        if texto == "" or texto.startswith("$"):
            continue

        campos = separar_campos(linea)

        if leyendo_nodos and len(campos) >= 4:
            try:
                nid = int(campos[0])
                x = float(campos[1])
                y = float(campos[2])
                z = float(campos[3])
                nodos[nid] = (x, y, z)

                if nid > max_nid:
                    max_nid = nid
            except:
                pass

        if leyendo_solidos and len(campos) >= 10:
            try:
                eid = int(campos[0])
                pid = int(campos[1])
                conectividad = [int(c) for c in campos[2:10]]
                elementos[eid] = {
                    "pid": pid,
                    "conectividad": conectividad,
                    "linea": i,
                }

                if eid > max_eid:
                    max_eid = eid
            except:
                pass

    return lineas, nodos, elementos, max_nid, max_eid


def existe_pid_en_part(lineas, pid_objetivo):
    for i, linea in enumerate(lineas):
        if not linea.strip().upper().startswith("*PART"):
            continue

        for j in range(i + 1, min(i + 8, len(lineas))):
            campos = separar_campos(lineas[j])

            if len(campos) >= 3:
                try:
                    pid = int(campos[0])

                    if pid == pid_objetivo:
                        return True
                except:
                    pass

    return False


def existe_mid_en_material(lineas, mid_objetivo):
    for i, linea in enumerate(lineas):
        if not linea.strip().upper().startswith("*MAT_"):
            continue

        for j in range(i + 1, min(i + 8, len(lineas))):
            campos = separar_campos(lineas[j])

            if len(campos) >= 2:
                try:
                    mid = int(campos[0])

                    if mid == mid_objetivo:
                        return True
                except:
                    pass

    return False


def obtener_mid_material(bloque):
    for linea in bloque[1:]:
        texto = linea.strip()

        if texto == "" or texto.startswith("$"):
            continue

        try:
            return int(linea[:10])
        except:
            pass

        campos = separar_campos(linea)

        if len(campos) == 0:
            continue

        try:
            return int(campos[0])
        except:
            pass

    return None


def reemplazar_materiales_celula(lineas):
    """
    Remueve bloques MAT existentes para MID 2 y MID 5 y agrega los nuevos.
    Esto permite parametrizar tanto citoplasma como nucleo.
    """
    nuevas = []
    i = 0
    removidos = []

    while i < len(lineas):
        if lineas[i].strip().upper().startswith("*MAT_"):
            j = i + 1

            while j < len(lineas) and not lineas[j].strip().startswith("*"):
                j += 1

            bloque = lineas[i:j]
            mid = obtener_mid_material(bloque)

            if mid in (mid_citoplasma, mid_nucleo):
                removidos.append(mid)
                i = j
                continue

        nuevas.append(lineas[i])
        i += 1

    bloque_materiales = []
    bloque_materiales.extend(
            mat_block(
            mid_citoplasma,
            ro_citoplasma,
            c1_citoplasma,
            xk_citoplasma,
            "Material_citoplasma",
        )
    )
    bloque_materiales.extend(
            mat_block(
            mid_nucleo,
            ro_nucleo,
            c1_nucleo,
            xk_nucleo,
            "Material_nucleo",
        )
    )

    nuevas = insertar_antes_de_keyword(nuevas, bloque_materiales, "*DEFINE_BOX_TITLE")

    return nuevas, removidos


def num_txt(valor):
    if valor == 0:
        return "0.0"

    if 1.0e-3 <= abs(valor) < 1.0e4:
        return f"{valor:.6g}"

    texto = f"{valor:.5E}"
    texto = texto.replace("E+0", "E").replace("E+", "E")
    texto = texto.replace("E-0", "E-")

    if len(texto) > 10:
        texto = f"{valor:.3E}"
        texto = texto.replace("E+0", "E").replace("E+", "E")
        texto = texto.replace("E-0", "E-")

    return texto


def linea_control_timestep(dtinit, tssfac, isdo, tslimt, dt2ms, lctm, erode, ms1st):
    return (
        f"{num_txt(dtinit):>10}"
        f"{num_txt(tssfac):>10}"
        f"{int(isdo):>10}"
        f"{num_txt(tslimt):>10}"
        f"{num_txt(dt2ms):>10}"
        f"{int(lctm):>10}"
        f"{int(erode):>10}"
        f"{int(ms1st):>10}\n"
    )


def linea_timestep_opt1():
    return (
        f"{0.0:10.1f}"
        f"{0:10d}"
        f"{1:10d}"
        f"{'':10}"
        f"{'':10}"
        f"{0.0:10.1f}"
        f"{0.0:10.1f}"
        f"{0:10d}\n"
    )


def leer_linea_timestep(linea):
    campos = separar_campos(linea)

    if len(campos) >= 8:
        return [
            float(campos[0]),
            float(campos[1]),
            int(float(campos[2])),
            float(campos[3]),
            float(campos[4]),
            int(float(campos[5])),
            int(float(campos[6])),
            int(float(campos[7])),
        ]

    return [0.0, 0.669, 0, 0.0, 0.0, 0, 0, 0]


def poner_dt2ms(lineas, dt2ms):
    comentario_opt1 = (
        "$#  dt2msf   dt2mslc     imscl    unused    unused"
        "     rmscl     emscl      ihdo\n"
    )

    for i, linea in enumerate(lineas):
        if not linea.strip().upper().startswith("*CONTROL_TIMESTEP"):
            continue

        j = i + 1

        while j < len(lineas):
            texto = lineas[j].strip()

            if texto == "" or texto.startswith("$"):
                j += 1
                continue

            if texto.startswith("*"):
                break

            vals = leer_linea_timestep(lineas[j])
            vals[4] = dt2ms
            lineas[j] = linea_control_timestep(*vals)

            k = j + 1

            while k < len(lineas) and not lineas[k].strip().startswith("*"):
                texto_opt = lineas[k].strip().upper()

                if "DT2MSF" in texto_opt or "IMSCL" in texto_opt:
                    lineas[k] = comentario_opt1
                    m = k + 1

                    while m < len(lineas) and lineas[m].strip() == "":
                        m += 1

                    if m < len(lineas) and not lineas[m].strip().startswith("*"):
                        lineas[m] = linea_timestep_opt1()
                    else:
                        lineas.insert(m, linea_timestep_opt1())

                    return True

                k += 1

            lineas[j + 1:j + 1] = [comentario_opt1, linea_timestep_opt1()]
            return True

    bloque = [
        "*CONTROL_TIMESTEP\n",
        "$#  dtinit    tssfac      isdo    tslimt     dt2ms      lctm     erode     ms1st\n",
        linea_control_timestep(0.0, 0.669, 0, 0.0, dt2ms, 0, 0, 0),
        comentario_opt1,
        linea_timestep_opt1(),
    ]
    nuevas = insertar_antes_de_keyword(lineas, bloque, "*DATABASE")
    lineas[:] = nuevas

    return True


# ============================================================
# 4. INTEGRACION
# ============================================================

lineas_modelo, nodos_modelo, elementos_modelo, max_nid, max_eid = leer_nodos_y_elementos_solidos(
    archivo_modelo_entrada
)

lineas_celula, nodos_celula, elementos_celula, max_nid_celula, max_eid_celula = leer_nodos_y_elementos_solidos(
    archivo_celula_nueva
)

mapa_nodos = {}
nodos_nuevos = []

nid_siguiente = max_nid + 1

for nid in sorted(nodos_celula):
    mapa_nodos[nid] = nid_siguiente
    nodos_nuevos.append((nid_siguiente, nodos_celula[nid]))
    nid_siguiente += 1

elementos_nuevos = []
eid_siguiente = max_eid + 1

for eid in sorted(elementos_celula):
    pid = elementos_celula[eid]["pid"]
    conectividad = elementos_celula[eid]["conectividad"]
    conectividad_nueva = [mapa_nodos[nid] for nid in conectividad]
    elementos_nuevos.append((eid_siguiente, pid, conectividad_nueva))
    eid_siguiente += 1

lineas_salida = []

leyendo_solidos = False
insertados_elementos_nuevos = False
elementos_pid2_removidos = 0

for linea in lineas_modelo:
    texto = linea.strip()

    if texto.upper().startswith("*ELEMENT_SOLID"):
        leyendo_solidos = True
        lineas_salida.append(linea)
        continue

    if leyendo_solidos and texto.startswith("*"):
        for eid, pid, conectividad in elementos_nuevos:
            lineas_salida.append(escribir_elemento(eid, pid, conectividad))

        insertados_elementos_nuevos = True
        leyendo_solidos = False
        lineas_salida.append(linea)
        continue

    if leyendo_solidos:
        if texto == "" or texto.startswith("$"):
            lineas_salida.append(linea)
            continue

        campos = separar_campos(linea)

        if len(campos) >= 10:
            try:
                pid = int(campos[1])

                if pid == pid_citoplasma:
                    elementos_pid2_removidos += 1
                    continue
            except:
                pass

    lineas_salida.append(linea)

if leyendo_solidos and not insertados_elementos_nuevos:
    for eid, pid, conectividad in elementos_nuevos:
        lineas_salida.append(escribir_elemento(eid, pid, conectividad))

if not insertados_elementos_nuevos:
    print("Advertencia: no se encontro cierre de *ELEMENT_SOLID; revisa el archivo.")


# Insertar nodos nuevos justo antes del primer *PART, ubicado despues del
# bloque *NODE en el modelo base. Asi no se modifica ningun nodo existente.
lineas_nodos_nuevos = ["$\n", "$ Nodos agregados para celula nucleo-citoplasma\n", "$\n"]

for nid, coordenadas in nodos_nuevos:
    lineas_nodos_nuevos.append(escribir_nodo(nid, coordenadas))

lineas_salida = insertar_antes_de_keyword(lineas_salida, lineas_nodos_nuevos, "*PART")

lineas_salida, materiales_removidos = reemplazar_materiales_celula(lineas_salida)

if not existe_pid_en_part(lineas_salida, pid_nucleo):
    lineas_salida = insertar_antes_de_keyword(
        lineas_salida,
        crear_bloque_part_nucleo(),
        "*CONTACT",
    )

if args.dt2ms is not None:
    poner_dt2ms(lineas_salida, args.dt2ms)


# ============================================================
# 5. GUARDAR Y REPORTAR
# ============================================================

archivo_modelo_salida.parent.mkdir(parents=True, exist_ok=True)

with open(archivo_modelo_salida, "w", encoding="utf-8") as f:
    f.writelines(lineas_salida)

print("Listo.")
print("Archivo generado:", archivo_modelo_salida)
print("")
print("Modelo base:", archivo_modelo_entrada)
print("Celula nueva:", archivo_celula_nueva)
print("")
print("Max node original:", max_nid)
print("Max elem original:", max_eid)
print("Nodos nuevos agregados:", len(nodos_nuevos))
print("Elementos solidos PID 2 antiguos removidos:", elementos_pid2_removidos)
print("Elementos nuevos agregados:", len(elementos_nuevos))
print("")
print("Materialidad:")
print("  preset:", materiales["preset"])
print("  materiales MID removidos/reemplazados:", materiales_removidos)
print("  Ec kPa:", materiales["Ec_kPa"])
print("  En kPa:", materiales["En_kPa"])
print("  En/Ec:", materiales["relacion_En_Ec"])
print("  nu citoplasma:", materiales["nu_citoplasma"])
print("  nu nucleo:", materiales["nu_nucleo"])
print("  C1 citoplasma:", c1_citoplasma)
print("  XK citoplasma:", xk_citoplasma)
print("  C1 nucleo:", c1_nucleo)
print("  XK nucleo:", xk_nucleo)
if args.dt2ms is not None:
    print("  DT2MS:", args.dt2ms)
print("")
print("PID citoplasma:", pid_citoplasma)
print("PID nucleo:", pid_nucleo)
