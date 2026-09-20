# Postproceso de campañas paramétricas

`analizar_campanas.py` extrae las coordenadas X de los nodos frontales del
citoplasma y del núcleo directamente desde los `d3plot`. Los nodos se
identifican en cada modelo `.k` como los de mayor coordenada X inicial de las
partes PID 2 (citoplasma) y PID 5 (núcleo).

La comparación usa:

```text
lambda_c = (x_c(t) - x_c(0)) / Rc
lambda_n = (x_n(t) - x_n(0)) / Rc
```

El mapa de calor informa `lambda_c` en un instante común. De forma
predeterminada se usa 0,035 ms, anterior al último estado legible más temprano de los
resultados actualmente disponibles.

## Uso

Desde `tesis_lsdyna`:

```powershell
python -m pip install -r scripts/postprocess/requirements.txt
python scripts/postprocess/analizar_campanas.py --campana todas
```

Para cambiar el instante de comparación:

```powershell
python scripts/postprocess/analizar_campanas.py --tiempo-referencia-ms 0.035
```

Las tablas y figuras se escriben en `resultados_procesados/`. Los casos sin
`d3plot` se conservan en las tablas con estado `pendiente` y aparecen como
celdas grises en los mapas; el script no interpola campañas ausentes.

## Salidas principales

- `material/resumen_material.csv`
- `material/trayectorias_material.csv`
- `material/figuras/material_trayectorias_normalizadas.{png,pdf}`
- `material/figuras/material_mapa_calor.{png,pdf}`
- `material/figuras/avance_normalizado_En100.{png,pdf}`
- `material/figuras/avance_normalizado_En225.{png,pdf}`
- `material/figuras/avance_normalizado_En300.{png,pdf}`
- `material/resumen_avance_constriccion.csv`
- `material/figuras/separacion_interna_En100.{png,pdf}`
- `material/figuras/separacion_interna_En225.{png,pdf}`
- `material/figuras/separacion_interna_En300.{png,pdf}`
- `material/resumen_separacion_interna.csv`
- `geometria/resumen_geometria.csv`
- `geometria/trayectorias_geometria.csv`
- `geometria/figuras/geometria_trayectorias_normalizadas.{png,pdf}`
- `geometria/figuras/geometria_mapa_calor.{png,pdf}`
- `geometria/figuras/progreso_constriccion_Rc0p0045.{png,pdf}`
- `geometria/figuras/progreso_constriccion_Rc0p0050.{png,pdf}`
- `geometria/figuras/progreso_constriccion_Rc0p0055.{png,pdf}`
- `geometria/resumen_progreso_constriccion.csv`

Los CSV de `series_lspp/` son la extracción cruda de LS-PrePost y permiten
auditar las tablas agregadas.

## Avance normalizado en la campaña de materialidad

`avance_materialidad.py` toma las trayectorias ya extraídas y genera una figura
independiente para cada valor de `En`:

```powershell
python scripts/postprocess/avance_materialidad.py
```

Para cada caso se representa
`xi_c = (x_c(t) - x_ent)/Rc` frente a `Delta t_ent = t - t_ent`, desde que el nodo
frontal del citoplasma ingresa a la constricción hasta que completa su
tránsito. Los instantes de cruce se interpolan linealmente entre los estados de
salida que encierran cada extremo. No se aplica suavizado ni filtrado a la
trayectoria.

## Separación interna núcleo-citoplasma

`separacion_materialidad.py` utiliza las mismas trayectorias para representar
el cambio de separación axial entre los nodos frontales del citoplasma y del
núcleo:

```powershell
python scripts/postprocess/separacion_materialidad.py
```

La magnitud graficada es
`Delta d_nc/Rc = [(x_c-x_n)-(x_c0-x_n0)]/Rc = lambda_c-lambda_n` y el eje
horizontal corresponde a `xi_c`. Las tres figuras comparten la misma escala
vertical. Los extremos del tránsito se interpolan linealmente y no se aplica
suavizado.

## Presión en la campaña de materialidad

`presion_materialidad.py` obtiene la presión de los elementos de agua
contenidos en dos láminas que encierran la constricción, calcula
`DeltaP = P_up - P_down` y representa la presión normalizada `DeltaP/Ec` frente
al tiempo transcurrido desde el ingreso a la constricción. Genera una figura
independiente por módulo nuclear:

```powershell
python scripts/postprocess/presion_materialidad.py
```

Las salidas se escriben en
`resultados_procesados/material/presion/`. Este análisis se conserva como
postproceso exploratorio y no se incorpora al capítulo de resultados mientras
se mantenga la condición de borde de flujo actualmente utilizada.

## Progreso en la constricción para la campaña geométrica

`progreso_geometria.py` representa el progreso del nodo frontal del citoplasma
dentro de la constricción:

```powershell
python scripts/postprocess/progreso_geometria.py
```

La magnitud vertical es
`chi = (x_c(t)-x_ent)/(x_sal-x_ent)` y el eje horizontal es
`Delta t_ent = t-t_ent`. Se genera una figura independiente por radio celular `Rc`,
con las tres relaciones geométricas `eta`. Los cruces de entrada y salida se
interpolan linealmente entre estados consecutivos y no se aplica suavizado.
Cuando un caso no alcanza la salida, la curva termina en el último estado
disponible, sin extrapolación.
