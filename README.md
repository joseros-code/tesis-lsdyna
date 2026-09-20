# Códigos de la tesis

Este repositorio reúne los scripts utilizados para construir y analizar el
modelo numérico de una célula multicomponente en un citómetro de constricción
mediante LS-DYNA.

## Contenido

```text
.
├── caso_base/
│   └── modelo_NC-08_En300_Ec150_Gamma2p00.k
└── scripts/
    ├── core/
    ├── campaigns/
    └── postprocess/
```

- `caso_base/` contiene el caso hiperelástico NC-08 empleado como referencia:
  $E_n=300$ kPa, $E_c=150$ kPa y $\Gamma_E=2$.
- `scripts/core/` contiene la definición de materiales, la generación de la
  malla celular y su integración en el modelo FSI-ALE.
- `scripts/campaigns/` contiene los barridos de sensibilidad de malla,
  materialidad, geometría y viscoelasticidad.
- `scripts/postprocess/` contiene las rutinas empleadas para obtener las
  métricas y figuras presentadas en la tesis.

## Requisitos

- Python 3.10 o posterior.
- LS-DYNA para ejecutar los modelos `.k`.
- LS-PrePost para la extracción automatizada de resultados desde `d3plot`.
- Las dependencias indicadas en `scripts/postprocess/requirements.txt`.

## Alcance

El repositorio no incluye los archivos binarios generados por LS-DYNA ni las
copias de todos los modelos paramétricos. Estos modelos se obtienen mediante
los scripts a partir del caso base y de los parámetros definidos en cada
barrido.

LS-DYNA y LS-PrePost son programas externos y no forman parte de este
repositorio.

## Autor

José Ignacio Ros Quezada  
Tesis de Ingeniería Civil Mecánica, Universidad de Santiago de Chile, 2026.

