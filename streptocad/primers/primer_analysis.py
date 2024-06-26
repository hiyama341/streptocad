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

import pandas as pd
from typing import List, Dict, Any
import primer3

def analyze_primers_and_hairpins(primer_df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyzes primer sequences for homodimers, heterodimers, and hairpins.

    Parameters
    ----------
    primer_df : pd.DataFrame
        DataFrame containing primer information with columns:
        - 'f_primer_sequences(5-3)'
        - 'r_primer_sequences(5-3)'
        - 'f_primer_name'
        - 'r_primer_name'
        - 'ta' (annealing temperature)

    Returns
    -------
    pd.DataFrame
        DataFrame containing analysis results for each primer pair with columns:
        - 'annealing_temperature'
        - 'primer_pair'
        - 'homodimer_forward_tm'
        - 'homodimer_forward_deltaG (kcal/mol)'
        - 'homodimer_reverse_tm'
        - 'homodimer_reverse_deltaG (kcal/mol)'
        - 'heterodimer_tm'
        - 'heterodimer_deltaG (kcal/mol)'
        - 'hairpin_forward_structure_found'
        - 'hairpin_forward_tm'
        - 'hairpin_forward_deltaG (kcal/mol)'
        - 'hairpin_reverse_structure_found'
        - 'hairpin_reverse_tm'
        - 'hairpin_reverse_deltaG (kcal/mol)'
    """
    results: List[Dict[str, Any]] = []

    for i, row in primer_df.iterrows():
        f_sequence: str = row['f_primer_sequences(5-3)']
        r_sequence: str = row['r_primer_sequences(5-3)']
        f_name: str = row['f_primer_name']
        r_name: str = row['r_primer_name']
        ta: float = row['ta']

        if len(f_sequence) > 60 or len(r_sequence) > 60:
            print(f"Skipping analysis for primer pair {f_name} & {r_name} due to length constraints.")
            continue

        try:
            homodimer_forward: Dict[str, Any] = primer3.bindings.calc_homodimer(f_sequence, temp_c=ta).todict()
            homodimer_reverse: Dict[str, Any] = primer3.bindings.calc_homodimer(r_sequence, temp_c=ta).todict()
            heterodimer: Dict[str, Any] = primer3.bindings.calc_heterodimer(f_sequence, r_sequence, temp_c=ta).todict()
            hairpin_forward: Dict[str, Any] = primer3.bindings.calc_hairpin(f_sequence, temp_c=ta).todict()
            hairpin_reverse: Dict[str, Any] = primer3.bindings.calc_hairpin(r_sequence, temp_c=ta).todict()
        except RuntimeError as e:
            print(f"Error during analysis for primer pair {f_name} & {r_name}: {e}")
            continue

        row_result: Dict[str, Any] = {
            "annealing_temperature": ta,
            "primer_pair": f"{f_name} & {r_name}",
            "homodimer_forward_tm": homodimer_forward['tm'],
            "homodimer_forward_deltaG (kcal/mol)": homodimer_forward['dg'] / 1000,
            "homodimer_reverse_tm": homodimer_reverse['tm'],
            "homodimer_reverse_deltaG (kcal/mol)": homodimer_reverse['dg'] / 1000,
            "heterodimer_tm": heterodimer['tm'],
            "heterodimer_deltaG (kcal/mol)": heterodimer['dg'] / 1000,
            "hairpin_forward_structure_found": hairpin_forward['structure_found'],
            "hairpin_forward_tm": hairpin_forward['tm'],
            "hairpin_forward_deltaG (kcal/mol)": hairpin_forward['dg'] / 1000,
            "hairpin_reverse_structure_found": hairpin_reverse['structure_found'],
            "hairpin_reverse_tm": hairpin_reverse['tm'],
            "hairpin_reverse_deltaG (kcal/mol)": hairpin_reverse['dg'] / 1000
        }
        results.append(row_result)

    df_results: pd.DataFrame = pd.DataFrame(results)
    return df_results




