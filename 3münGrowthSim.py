#!/usr/bin/env python
# coding: utf-8

# In[2]:


#!python3 
"""

@authors: Felipe Ulloa-Fierro & Fernando Badilla Véliz
v2.0

"""

import copy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import List

# Cargar el archivo CSV
df = pd.read_csv("lookup_table.csv", keep_default_na=False)
df.set_index(
    ["Especie", "Zona", "SiteIndex", "Manejo", "Condicion", "DensidadInicial"],
    inplace=True,
)
dict_idx = df["id"].to_dict()  # growth_equations_id_dictionary

config = {
    "horizonte": 40,  # Provisto por el usuario
    "rodales": 4,  # Proviene del shapefile
    "edades": [5, 4, 6, 10],  # Proviene del shapefile
    "has": [2, 3, 4, 5],  # Proviene del shapefile
    "especie": ["Pinus", "Pinus", "Pinus", "Eucapyltus"],  # Proviene del shapefile
    "siteIndex": [32, 29, 23, 30],  # Proviene del shapefile
    "zona": [6, 6, 6, 2],  # Proviene del shapefile
    "densidadInicial": [1250, 1250, 1250, 1250],
    "manejo": ["Intensivo", "Intensivo2", "Pulpable", "NA"],  # Proviene del shapefile
    "condición": [
        "PostRaleo1250-700",
        "PostRaleo1250-700",
        "PostRaleo1250-700",
        "SinManejo",
    ],  # Proviene del shapefile
    "id_rodal": ["A11", "D22", "C33", "F34"],  # Proviene del shapefile
    "num_policies_pino": 3,  # Provisto por el usuario
    "policies_pino": [[8, 25], [10, 27], [12, 30]],  # Formato: [Periodo Poda_Raleo, Periodo Cosecha]  # Provisto por el usuario
    "num_policies_eucalyptus": 2,  # Provisto por el usuario
    "policies_eucalyto": [[8], [10]],  # años de cosecha ya que en eucalyptus no se ralea.  # Provisto por el usuario
}


def calc_biomasa(rodal: pd.Series, e: int) -> float:
    return max(rodal["α"] * e ** rodal["β"] + rodal["γ"], 0)


def generar_rodales(
    edades=config["edades"],
    has=config["has"],
    id_rodal=config["id_rodal"],
    tipoespecie=config["especie"],
) -> List[pd.Series]:
    """Generar una lista de rodales usando la configuración dada"""
    rodales = []
    if "id" in df.columns:
        df.set_index("id", inplace=True)  # para asegurarse que 'id' es el índice de df
    for i in range(config["rodales"]):
        key = (
            config["especie"][i],
            config["zona"][i],
            config["siteIndex"][i],
            config["manejo"][i],
            config["condición"][i],
            config["densidadInicial"][i],
        )
        if key not in dict_idx:
            raise KeyError(
                f"No existe un índice para la combinación de claves {key} en la lookup table."
            )
        idx = dict_idx[key]
        if idx not in df.index:  # checkea si idx está en df.index en lugar de df['id'].values
            raise ValueError(f"El índice {idx} no existe en la lookup table.")
        init_ages = edades[i]
        hectareas = has[i]
        rodal_id = id_rodal[i]
        especie = tipoespecie[i]
        rodal = pd.concat(
            (
                df.loc[idx],
                pd.Series(
                    {
                        "eq_id": idx,
                        "edad_in": init_ages,
                        "ha": hectareas,
                        "id_rodal": rodal_id,
                        "TipoEspecie": especie,
                    }
                ),
            )
        )
        rodales.append(rodal)
    return rodales


def generar_rodalesconpolicy(
    rodales=generar_rodales(), horizonte=config["horizonte"]
):
    """Generar una lista de rodales con políticas de manejo aplicadas"""
    rodales_con_policy = []
    for rodal in rodales:
        if rodal["TipoEspecie"] == "Pinus":
            policies = config["policies_pino"]
        elif rodal["TipoEspecie"] == "Eucapyltus":
            policies = config["policies_eucalyto"]
        else:
            raise ValueError(
                f"La especie {rodal['TipoEspecie']} no tiene políticas definidas."
            )

        for policy in policies:
            rodal_policy = rodal.copy()
            poda_raleo_periodos = []
            cosecha_periodos = []
            periodo_actual = 0
            while periodo_actual <= horizonte:
                if rodal["TipoEspecie"] == "Pinus":
                    periodo_poda_raleo = periodo_actual + policy[0]
                    periodo_cosecha = periodo_actual + policy[1]
                    if periodo_poda_raleo <= horizonte:
                        poda_raleo_periodos.append(periodo_poda_raleo)
                    if periodo_cosecha <= horizonte:
                        cosecha_periodos.append(periodo_cosecha)
                    periodo_actual = periodo_cosecha
                elif rodal["TipoEspecie"] == "Eucapyltus":
                    periodo_cosecha = periodo_actual + policy[0]
                    if periodo_cosecha <= horizonte:
                        cosecha_periodos.append(periodo_cosecha)
                    periodo_actual = periodo_cosecha
            rodal_policy["poda_raleo"] = poda_raleo_periodos
            rodal_policy["cosecha"] = cosecha_periodos
            rodal_policy["horizonte"] = horizonte

            rodales_con_policy.append(rodal_policy)

    return rodales_con_policy



def rotular_condicion(
    edad_poda_raleo: List[int], edad_cosecha: List[int], edades_updated: pd.Series
) -> pd.Series:
    condicion = ["sin manejo"] * len(edades_updated)

    if len(edad_poda_raleo) == 0 and len(edad_cosecha) == 0:
        return pd.Series(condicion)

    indice_poda_raleo = 0
    indice_cosecha = 0
    en_manejo = False

    for i in range(len(edades_updated)):
        edad = edades_updated[i]

        if indice_poda_raleo < len(edad_poda_raleo) and edad == edad_poda_raleo[indice_poda_raleo]:
            en_manejo = True
            indice_poda_raleo += 1

        if indice_cosecha < len(edad_cosecha) and edad == edad_cosecha[indice_cosecha]:
            en_manejo = False
            indice_cosecha += 1

        if en_manejo:
            condicion[i] = "con manejo"

    return pd.Series(condicion)


def generar_codigo_kitral(especie: str, edad: int, condición: str) -> str:
    """Genera un diccionario de códigos Kitral basado en la especie, edad y condición"""
    key = (especie, edad, condición)
    if key[0] == "Pinus":
        if key[1] <= 3:
            value = "PL01"
        elif 3 < key[1] <= 11:
            value = "PL05" if key[2] == "con manejo" else "PL02"
        elif 11 < key[1] <= 17:
            value = "PL06" if key[2] == "con manejo" else "PL03"
        else:
            value = "PL07" if key[2] == "con manejo" else "PL04"
    else:  # Eucalyptus
        if key[1] <= 3:
            value = "PL08"
        elif 3 < key[1] <= 10:
            value = "PL09"
        else:
            value = "PL10"
    return value

def simula_bosque(
    rodales_con_policy=generar_rodalesconpolicy(),
    horizonte=config["horizonte"],
    num_policies_pino=config["num_policies_pino"],
    num_policies_eucalyptus=config["num_policies_eucalyptus"],
):
    """Simula el crecimiento del bosque calculando la biomasa para cada rodal con política"""
    bosque = []
    cantidad_raleada_total = []
    cantidad_cosechada_total = []
    resumen = []
    # Inicializa los contadores de políticas
    policy_counters = {
        "Pinus": 0,
        "Eucapyltus": 0,
    }

    for rodal in rodales_con_policy:
        tabla = pd.DataFrame()
        tabla["periodo"] = range(1, horizonte + 1)
        tabla["edadRodal"] = range(rodal["edad_in"], rodal["edad_in"] + horizonte)

        # para crear la columna "edad_updated"
        tabla["edad_updated"] = tabla["edadRodal"]

        # encontramos los índices donde ocurren las cosechas
        cosecha_indices = [
            tabla[tabla["periodo"] == c].index[0] for c in rodal["cosecha"]
        ]

        # refresh de los valores de "edad_updated" después de cada cosecha
        for cosecha_index in cosecha_indices:
            longitud = len(
                tabla.loc[cosecha_index + 1 :, "edad_updated"]
            )  # comienza desde el periodo siguiente a la cosecha
            valores = range(1, longitud + 1)
            tabla.loc[
                cosecha_index + 1 :, "edad_updated"
            ] = valores  # actualiza desde el periodo siguiente a la cosecha

        # define edad_poda_raleo y edad_cosecha aquí
        edad_poda_raleo = [
            tabla.loc[tabla["periodo"] == pr, "edad_updated"].values[0]
            for pr in rodal["poda_raleo"]
        ] if rodal["poda_raleo"] is not None else []
        edad_cosecha = [
            tabla.loc[tabla["periodo"] == pc, "edad_updated"].values[0]
            for pc in rodal["cosecha"]
        ]

        # Rotular la columna "condición"
        tabla["condición"] = rotular_condicion(
            edad_poda_raleo, edad_cosecha, tabla["edad_updated"]
        )

        if pd.isnull(rodal.get("next")) or rodal["next"] == "":
            rodal["poda_raleo"] = -1
            biomasa = []
            for e in tabla["edad_updated"]:
                biomasa.append(calc_biomasa(rodal, e))
        else:
            next_rodal = df.loc[int(rodal["next"])]
            biomasa = []
            poda_raleo_index = 0
            cosecha_index = 0

            for e in tabla["edad_updated"]:
                if (
                    poda_raleo_index < len(edad_poda_raleo)
                    and e == edad_poda_raleo[poda_raleo_index]
                ):
                    biomasa.append(calc_biomasa(rodal, e))
                    poda_raleo_index += 1
                elif (
                    cosecha_index < len(edad_cosecha)
                    and e == edad_cosecha[cosecha_index]
                ):
                    biomasa.append(calc_biomasa(next_rodal, e))
                    cosecha_index += 1
                else:
                    if (
                        poda_raleo_index > 0
                        and e > edad_poda_raleo[poda_raleo_index - 1]
                    ):
                        biomasa.append(calc_biomasa(next_rodal, e))
                    else:
                        biomasa.append(calc_biomasa(rodal, e))

        tabla["biomasa"] = (pd.Series(biomasa) * rodal["ha"]).round(3)
        tabla["id_rodal"] = rodal["id_rodal"]
        tabla["Especie"] = rodal["TipoEspecie"]

        # Crea la columna "kitral_class" utilizando la función generar_codigo_kitral
        tabla["kitral_class"] = tabla.apply(
            lambda row: generar_codigo_kitral(
                row["Especie"], row["edad_updated"], row["condición"]
            ),
            axis=1,
        )

        edad_initial = tabla.loc[
            tabla["edadRodal"] == rodal.edad_in, "edad_updated"
        ].values[0]
        edad_final = tabla.loc[
            tabla["periodo"] == rodal.horizonte, "edad_updated"
        ].values[0]

        # agrega los valores a la tabla de resumen
        if rodal["TipoEspecie"] == "Pinus":
            politica = f"policy_pino {policy_counters['Pinus'] % num_policies_pino + 1}"
            policy_counters["Pinus"] += 1
        elif rodal["TipoEspecie"] == "Eucapyltus":
            politica = f"policy_eucalyptus {policy_counters['Eucapyltus'] % num_policies_eucalyptus + 1}"
            policy_counters["Eucapyltus"] += 1
        else:
            raise ValueError(
                f"La especie {rodal['TipoEspecie']} no tiene políticas definidas."
            )

        tabla["politica"] = politica

        # Calcula la cantidad cosechada por cada tiempo de cosecha por rodal y policy
        cantidad_cosechada = []
        for pc in rodal["cosecha"]:
            if len(tabla.loc[tabla["periodo"] == pc, "biomasa"].values) > 0:
                cantidad_cosechada.append(
                    tabla.loc[tabla["periodo"] == pc, "biomasa"].values[0]
                )
            else:
                cantidad_cosechada.append(0)

        cantidad_cosechada_total.append(cantidad_cosechada)

        cantidad_raleada = []
        if isinstance(rodal["poda_raleo"], list):
            for pr in rodal["poda_raleo"]:
                if pr != -1:
                    if (
                        calc_biomasa(
                            next_rodal,
                            tabla.loc[tabla["periodo"] == pr, "edad_updated"].values[0],
                        )
                        > 0
                    ):
                        cantidad_raleada.append(
                            calc_biomasa(
                                rodal,
                                tabla.loc[tabla["periodo"] == pr, "edad_updated"].values[0],
                            )
                            - calc_biomasa(
                                next_rodal,
                                tabla.loc[tabla["periodo"] == pr, "edad_updated"].values[0],
                            )
                        )
                    else:
                        cantidad_raleada.append(
                            (
                                tabla.loc[tabla["periodo"] == pr, "biomasa"].values[0]
                                - tabla.loc[tabla["periodo"] == pr + 1, "biomasa"].values[0]
                            )
                        )
                else:
                    cantidad_raleada.append(0)

        cantidad_raleada_total.append(cantidad_raleada)

        cantidad_cosechada_total = [
            [round(x, 3) for x in sublist] for sublist in cantidad_cosechada_total
        ]
        cantidad_raleada_total = [
            [round(x, 3) for x in sublist] for sublist in cantidad_raleada_total
        ]

        bosque.append(tabla)

        resumen.append(
            {
                "id_rodal": rodal["id_rodal"],
                "especie": rodal["TipoEspecie"],
                "has": rodal["ha"],
                "edad_in": rodal["edad_in"],
                "edadin_kitralcode": tabla.loc[
                    tabla["edad_updated"] == edad_initial, "kitral_class"
                ].values[0],
                "politica": politica,
                "edad_poda_raleo": edad_poda_raleo,
                "cantidad(es)_raleada": [round(x, 3) for x in cantidad_raleada],
                "edad_cosecha": edad_cosecha,
                "cantidad(es)_cosechada": [round(x, 3) for x in cantidad_cosechada],
                "edad_final": edad_final,
                "edadfin_kitralcode": tabla.loc[
                    tabla["edad_updated"] == edad_final, "kitral_class"
                ].values[0],
            }
        )

    # Convierte la lista de resumen en un DataFrame
    resumen = pd.DataFrame(resumen)

    return bosque, resumen



bosque, resumen = simula_bosque(
    rodales_con_policy=generar_rodalesconpolicy(),
    horizonte=config["horizonte"],
    num_policies_pino=config["num_policies_pino"],
    num_policies_eucalyptus=config["num_policies_eucalyptus"],
)

# Exporta el DataFrame resumen a un archivo Excel
resumen.to_excel("forest_summary.xlsx", index=False)

# Separa el DataFrame bosque en dos DataFrames, uno para cada especie
bosque_pino = pd.concat([df for df in bosque if df["Especie"].iloc[0] == "Pinus"])
bosque_eucalyptus = pd.concat(
    [df for df in bosque if df["Especie"].iloc[0] == "Eucapyltus"]
)

# Para cada especie, crea un archivo Excel con una pestaña para cada política
with pd.ExcelWriter("bosque_pino.xlsx") as writer:
    for i in range(1, config["num_policies_pino"] + 1):
        df = bosque_pino[bosque_pino["politica"] == f"policy_pino {i}"]
        df.to_excel(writer, sheet_name=f"policy_pino {i}", index=False)

with pd.ExcelWriter("bosque_eucalyptus.xlsx") as writer:
    for i in range(1, config["num_policies_eucalyptus"] + 1):
        df = bosque_eucalyptus[
            bosque_eucalyptus["politica"] == f"policy_eucalyptus {i}"
        ]
        df.to_excel(writer, sheet_name=f"policy_eucalyptus {i}", index=False)


# In[ ]:




