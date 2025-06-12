# #!/usr/bin/env python
# # MIT License
# # Copyright (c) 2024, Technical University of Denmark (DTU)
# #
# # Permission is hereby granted, free of charge, to any person obtaining a copy
# # of this software and associated documentation files (the "Software"), to deal
# # in the Software without restriction, including without limitation the rights
# # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# # copies of the Software, and to permit persons to whom the Software is
# # furnished to do so, subject to the following conditions:
# #
# # The above copyright notice and this permission notice shall be included in all
# # copies or substantial portions of the Software.

import pandas as pd
from Bio.Seq import Seq


def identify_base_editing_sites(
    sgrna_df: pd.DataFrame, editing_window_start: int = 3, editing_window_end: int = 10
) -> pd.DataFrame:
    """
    Identify potential base editing sites for C-to-T substitutions in sgRNAs.
    """

    def find_editable_bases(
        row,
        editing_window_start=editing_window_start,
        editing_window_end=editing_window_end,
    ):
        sgrna = row["sgrna"]
        positions = []
        for i in range(editing_window_start, editing_window_end):
            if i < len(sgrna) and sgrna[i] == "C":
                positions.append(i + 1)

        return ",".join(map(str, positions))

    def find_context_dependent_seqs(row):
        sgrna = row["sgrna"]
        # gene_strand = row["gene_strand"]
        sequence_context_bases = []
        for i in range(editing_window_start, editing_window_end):
            if i < len(sgrna):
                # if gene_strand == 1:
                if sgrna[i] == "C":
                    # context: G immediately 5' of C
                    if i > 0 and sgrna[i - 1] == "G":
                        sequence_context_bases.append(1)
                    else:
                        sequence_context_bases.append(0)
            # elif gene_strand == -1:
            #     if sgrna[i] == "G":
            #         # context: C immediately 5' of G
            #         if i > 0 and sgrna[i - 1] == "C":
            #             sequence_context_bases.append(1)
            #         else:
            #             sequence_context_bases.append(0)
        return ",".join(map(str, sequence_context_bases))

    sgrna_df = sgrna_df.copy()
    sgrna_df["editing_context"] = sgrna_df.apply(find_context_dependent_seqs, axis=1)
    sgrna_df["editable_cytosines"] = sgrna_df.apply(find_editable_bases, axis=1)
    return sgrna_df


def filter_sgrnas_for_base_editing(sgrna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter sgRNAs to include only those with editable cytosines within the editing window.

    Parameters
    ----------
    sgrna_df : pd.DataFrame
        DataFrame containing sgRNA information.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with sgRNAs suitable for base editing.

    Examples
    --------
    >>> data = {'sgrna': ['GACCGT', 'CCGTGA'], 'editable_cytosines': ['3', '']}
    >>> df = pd.DataFrame(data)
    >>> result = filter_sgrnas_for_base_editing(df)
    >>> print(result)
    """
    return sgrna_df[sgrna_df["editable_cytosines"] != ""]


def process_base_editing(
    df: pd.DataFrame,
    gene_sequences: dict,
    only_stop_codons: bool = False,
    editing_context: bool = True,
) -> pd.DataFrame:
    """
    Process the DataFrame to apply C-to-T (or C-to-A if sgrna_strand is -1) mutations at specified positions and identify amino acid changes.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing sgRNA information with editable cytosines.
    gene_sequences : dict
        Dictionary mapping locus tags to gene sequences.
    only_stop_codons : bool
        If True, filter out rows without STOP codon mutations.


    Returns
    -------
    pd.DataFrame
        Updated DataFrame with new columns for mutated sequences and amino acid changes.
    """

    def mutate_sequence(row, gene_sequences):
        locus = row["locus_tag"]
        gene_seq = str(gene_sequences[locus])  # this is always 5'→3' of the CDS
        sgrna = row["sgrna"]  # your 5'→3' protospacer
        editable = list(map(int, row["editable_cytosines"].split(",")))
        strand = row["sgrna_strand"]  # +1 or -1, from your DataFrame

        # 1) pick the orientation of the protospacer so it matches gene_seq
        if strand == row["gene_strand"]:
            protospacer = sgrna
        else:
            protospacer = str(Seq(sgrna).reverse_complement())

        # 2) find where that 20mer actually sits in the gene
        start_idx = gene_seq.find(protospacer)
        if start_idx < 0:
            raise RuntimeError(f"Can't find {protospacer} in gene {locus}")

        # turn gene_seq into a mutable list
        mutated = list(gene_seq)

        # 3) for each editable base convert at the correct offset
        for pos in editable:
            if strand == row["gene_strand"]:
                # left‐to‐right mapping for +/+
                gpos = start_idx + pos - 1
            else:
                # if you reversed the sgRNA to match gene_seq,
                # the first sgRNA base sits at start_idx, so
                # a base at “pos” from the left of the original
                # sgRNA is at start_idx + (len- pos) in gene_seq

                gpos = start_idx + (len(sgrna) - pos)

            # sanity check bounds
            if not (0 <= gpos < len(mutated)):
                print(f"warning: out‐of‐bounds {gpos} in {locus}")
                continue

            # perform the edit
            ref = mutated[gpos]
            if ref == "C":
                mutated[gpos] = "T"
            elif ref == "G":
                mutated[gpos] = "A"
            else:
                # nothing to do
                continue

        # # build your string (e.g. from your `mutated` list)
        # nuc_seq = "".join(mutated)
        # # translate (default table, stops=*)
        # prot = Seq(nuc_seq).translate()
        # if "*" in prot[:-2]:
        #     print(
        #         f"{locus} {ref}{gpos}→{mutated[gpos]} {protospacer} {row['sgrna_strand']}"
        #     )
        #     print("".join(gene_seq))
        #     print(nuc_seq)

        # print(prot)
        # # print(
        # #     f"{locus} {ref}{gpos}→{mutated[gpos]} {protospacer} {row['sgrna_strand']}"
        # # )
        # # # print the original and mutated sequences for debugging
        # # print("".join(list(gene_seq)))
        # # print("".join(mutated))

        return "".join(mutated)

    def translate_sequence(sequence):
        """
        Translate a nucleotide sequence into an amino acid sequence.
        """
        return str(Seq(sequence).translate())

    def find_amino_acid_changes(row):
        original_seq = gene_sequences[row["locus_tag"]]
        mutated_seq = row["mutated_sequence"]

        original_aa_seq = translate_sequence(original_seq)
        mutated_aa_seq = translate_sequence(mutated_seq)

        mutations = []
        for i, (orig_aa, mut_aa) in enumerate(
            zip(original_aa_seq, mutated_aa_seq), start=0
        ):
            if orig_aa != mut_aa:
                mutations.append(f"{orig_aa}{i + 1}{mut_aa}")

        return ", ".join(mutations)

    def extract_first_mutation_position(mutations):
        if mutations:
            # Extract the position of the first mutation
            first_mutation = mutations.split(",")[0]
            position = int("".join(filter(str.isdigit, first_mutation)))
            return position
        return float("inf")

    # Use .loc to avoid SettingWithCopyWarning
    df = df.copy()
    df["mutated_sequence"] = df.apply(
        lambda row: mutate_sequence(row, gene_sequences), axis=1
    )
    df.loc[:, "mutations"] = df.apply(find_amino_acid_changes, axis=1)

    # Remove rows where mutations column is empty
    df = df[df["mutations"] != ""]
    # Drop the mutated_sequence column if you don't want it in the final output
    df = df.drop(columns=["mutated_sequence"])

    if only_stop_codons:
        # Filter out rows where there is no stop codon in the mutations column
        df = df[df["mutations"].str.contains(r"\*")]

        # Add a column for the position of the first mutation
        df["first_mutation_position"] = df["mutations"].apply(
            extract_first_mutation_position
        )

        # Sort the DataFrame by the position of the first mutation in ascending order
        df = df.sort_values(by="first_mutation_position", ascending=True)

        # Drop the first_mutation_position column if you don't want it in the final output
        df = df.drop(columns=["first_mutation_position"])

    if editing_context:
        # Filter out rows where 'editing_context' contains '1'
        df = df[~df["editing_context"].str.contains(r"1")]

    return df
