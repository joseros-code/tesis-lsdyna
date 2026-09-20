# ------------------------------------------------------------
# Materiales de la celula para *MAT_SOFT_TISSUE.
#
# Unidades consistentes del modelo base:
#   longitud: mm
#   masa: tonelada
#   tiempo: s
#   esfuerzo: MPa
#
# Por comodidad, los modulos de Young se ingresan en kPa y se convierten
# internamente a MPa.
# ------------------------------------------------------------


def mat_desde_E(E_kPa, nu):
    """
    Convierte E y nu a C1 y K. E entra en kPa y sale en MPa.
    """
    if nu >= 0.5:
        raise ValueError("nu debe ser menor que 0.5")

    if E_kPa <= 0.0:
        raise ValueError("E_kPa debe ser positivo")

    E_MPa = E_kPa / 1000.0
    c1 = E_MPa / (4.0 * (1.0 + nu))
    xk = E_MPa / (3.0 * (1.0 - 2.0 * nu))

    return c1, xk


def leer_mat(args):
    """
    Lee los datos de material desde los argumentos del script.
    """
    preset = args.material_preset

    if preset == "modelo_base":
        Ec_kPa = args.Ec_kPa if args.Ec_kPa is not None else 300.0
        relacion = args.relacion_En_Ec if args.relacion_En_Ec is not None else 1.0

    elif preset == "paper_tcell":
        Ec_kPa = args.Ec_kPa if args.Ec_kPa is not None else 0.078
        relacion = args.relacion_En_Ec if args.relacion_En_Ec is not None else 30.0

    elif preset == "manual":
        Ec_kPa = args.Ec_kPa if args.Ec_kPa is not None else 300.0
        relacion = args.relacion_En_Ec if args.relacion_En_Ec is not None else 1.0

    else:
        raise ValueError(f"Preset de material no reconocido: {preset}")

    En_kPa = args.En_kPa

    if En_kPa is None:
        En_kPa = Ec_kPa * relacion
    else:
        relacion = En_kPa / Ec_kPa

    c1_citoplasma, xk_citoplasma = mat_desde_E(
        Ec_kPa,
        args.nu_citoplasma,
    )
    c1_nucleo, xk_nucleo = mat_desde_E(
        En_kPa,
        args.nu_nucleo,
    )

    return {
        "preset": preset,
        "ro_citoplasma": args.ro_citoplasma,
        "ro_nucleo": args.ro_nucleo,
        "Ec_kPa": Ec_kPa,
        "En_kPa": En_kPa,
        "relacion_En_Ec": relacion,
        "nu_citoplasma": args.nu_citoplasma,
        "nu_nucleo": args.nu_nucleo,
        "c1_citoplasma": c1_citoplasma,
        "xk_citoplasma": xk_citoplasma,
        "c1_nucleo": c1_nucleo,
        "xk_nucleo": xk_nucleo,
    }


def args_mat(parser):
    parser.add_argument(
        "--material-preset",
        "--mat",
        choices=["modelo_base", "paper_tcell", "manual"],
        dest="material_preset",
        default="modelo_base",
        help="Preset de materialidad para citoplasma y nucleo.",
    )
    parser.add_argument(
        "--Ec-kPa",
        "--Ec",
        dest="Ec_kPa",
        type=float,
        default=None,
        help="Modulo de Young del citoplasma/resto celular en kPa.",
    )
    parser.add_argument(
        "--En-kPa",
        "--En",
        dest="En_kPa",
        type=float,
        default=None,
        help="Modulo de Young del nucleo en kPa. Si se omite usa En/Ec.",
    )
    parser.add_argument(
        "--relacion-En-Ec",
        "--EnEc",
        dest="relacion_En_Ec",
        type=float,
        default=None,
        help="Relacion En/Ec. Se ignora si se entrega --En-kPa.",
    )
    parser.add_argument(
        "--nu-citoplasma",
        "--nu-c",
        type=float,
        dest="nu_citoplasma",
        default=0.49999,
        help="Modulo de Poisson del citoplasma.",
    )
    parser.add_argument(
        "--nu-nucleo",
        "--nu-n",
        type=float,
        dest="nu_nucleo",
        default=0.49999,
        help="Modulo de Poisson del nucleo.",
    )
    parser.add_argument(
        "--ro-citoplasma",
        "--ro-c",
        type=float,
        dest="ro_citoplasma",
        default=1.4e-9,
        help="Densidad del citoplasma en unidades consistentes.",
    )
    parser.add_argument(
        "--ro-nucleo",
        "--ro-n",
        type=float,
        dest="ro_nucleo",
        default=1.4e-9,
        help="Densidad del nucleo en unidades consistentes.",
    )


def mat_block(mid, ro, c1, xk, titulo):
    return [
        "$\n",
        f"$ {titulo}\n",
        "$\n",
        "*MAT_SOFT_TISSUE_TITLE\n",
        f"{titulo}\n",
        "$#     mid        ro        c1        c2        c3        c4        c5\n",
        f"{mid:10d}{ro:10.4E}{c1:10.4E}"
        f"{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}\n",
        "$#      xk      xlam      fang     xlam0    failsf    failsm   failshr\n",
        f"{xk:10.4E}{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}"
        f"{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}\n",
        "$#    aopt        ax        ay        az        bx        by        bz\n",
        f"{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}"
        f"{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}\n",
        "$#     la1       la2       la3      macf\n",
        f"{0.0:10.4E}{0.0:10.4E}{0.0:10.4E}{1:10d}\n",
    ]

parametros_soft_tissue_desde_E = mat_desde_E
resolver_materiales = leer_mat
agregar_argumentos_material = args_mat
bloque_mat_soft_tissue = mat_block
