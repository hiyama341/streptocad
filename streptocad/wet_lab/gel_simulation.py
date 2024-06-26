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

from typing import List, Any
from pydna.gel import gel
from pydna.ladders import GeneRuler_1kb_plus


def simulate_gel_electrophoresis(list_of_amplicons: List[Any]) -> Any:
    """
    Simulates gel electrophoresis for a given list of amplicons.

    Parameters
    ----------
    list_of_amplicons : List[Any]
        A list of amplicon objects to be analyzed in the gel electrophoresis simulation.

    Returns
    -------
    Any
        A simulated gel electrophoresis image showing the migration of the amplicons through an agarose gel.

    Notes
    -----
    The function prints the name of each amplicon and then simulates a gel electrophoresis
    using these amplicons along with the GeneRuler 1kb Plus DNA ladder. The simulation
    visualizes how the amplicons would migrate through an agarose gel under an electric field.

    Examples
    --------
    >>> from pydna.amplify import Anneal
    >>> # Assuming `primer1`, `primer2`, and `template_dna` are defined pydna Primer and Dseqrecord objects
    >>> amplicon = Anneal(primers=[primer1, primer2], template=template_dna).products[0]
    >>> simulate_gel_electrophoresis([amplicon])
    """
    bands = list_of_amplicons
    for band in bands:
        print(band.name)

    return gel([GeneRuler_1kb_plus, *[[band] for band in bands]])
