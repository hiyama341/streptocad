#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

from teemi.build.PCR import (
    calculate_volumes,
)
from teemi.build.transformation import wanted_mass, wanted_volume
import pandas as pd

def calculate_master_mix(
    vol_p_reac=0, no_of_reactions=1, standard_reagents=[], standard_volumes=[]
):
    """Calculates the mastermix for number of reactions excluding primers"""

    # calculate the mix
    master_mix = calculate_volumes(
        vol_p_reac=vol_p_reac,
        no_of_reactions=no_of_reactions,
        standard_reagents=standard_reagents,
        standard_volumes=standard_volumes,
    )

    # Ensure the column that will hold mixed types is of object dtype
    master_mix.columns = ["vol_p_reac", f"mastermix for {no_of_reactions} reactions"]
    master_mix[f"mastermix for {no_of_reactions} reactions"] = master_mix[f"mastermix for {no_of_reactions} reactions"].astype(object)
    
    for k, v in master_mix.iterrows():
        if "Primer" in k or "primer" in k:
            master_mix.loc[
                k, f"mastermix for {no_of_reactions} reactions"
            ] = "Add primers individually"

    return master_mix

def calculate_volume_and_total_concentration_df(
    amplicons, amplicon_parts_amounts_total, n=1
):
    """Calculates the volume and total concentration of
    a list of DNA parts.

    Parameters
    ----------
    amplicons : list
        A list of amplicon objects
    amplicon_parts_amounts_total : dict
        A dictionary of amplicon names and their respective total amounts
    n : int (optional)
        Gives the option of multiplying the volume is needed. Optional set to 1.

    Returns
    -------
    volumes : list
        List of volumes of each amplicon
    ngs : list
        List of ngs of each amplicon
    total_conc : float
        Total concentration of all amplicons
    """
    volume_and_concentraion_list = []

    volumes = []
    ngs = []
    for amp in amplicons:
        w_moles = amplicon_parts_amounts_total[amp.name]
        w_mass = wanted_mass(wanted_moles=w_moles, size=len(amp))
        act_conc = amp.annotations["batches"][0]["concentration"]
        w_volume = wanted_volume(w_mass, act_conc) * n
        volumes.append(w_volume)
        ngs.append(w_volume * act_conc)

        volume_and_concentraion = {
            "name": amp.name,
            "volume to add": w_volume,
            "concentration": act_conc,
            "location": amp.annotations["batches"][0]["location"],
        }
        volume_and_concentraion_list.append(volume_and_concentraion)

    total_vol = sum(volumes)
    total_ngs = sum(ngs)
    total_conc = total_ngs / total_vol
    df = pd.DataFrame(volume_and_concentraion_list)
    print(f"Total volume of the parts mixed : {sum(volumes)}")
    print(f"Final concentration of the parts mixed : {total_conc}")

    return df


def dilute_solution(
    target_volume: float, concentration: float, final_concentration: float
):
    total_nanograms = target_volume * final_concentration

    solute_to_add = total_nanograms / concentration
    water_to_add = target_volume - solute_to_add

    return round(solute_to_add, 2), round(water_to_add, 2)


def calculate_BsaI_volume(number_of_fragments: int):
    BsaI_volume = 0.2 * number_of_fragments
    if BsaI_volume > 1:
        raise ValueError(
            f"In total, the amount of FastDigest BsaI should not exceed 1 µL, and you entered {BsaI_volume}"
        )

    return 0.2 * number_of_fragments
