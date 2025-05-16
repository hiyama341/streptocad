import pandas as pd
from typing import List
from Bio.SeqRecord import SeqRecord
from teemi.build.PCR import primer_tm_neb, primer_ta_neb
from pydna.design import primer_design
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.dseqrecord import Dseqrecord
import warnings
import pandas as pd
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design

from ..primers.primer_analysis import analyze_primers_and_hairpins
from ..cloning.pcr_simulation import make_amplicons


def generate_primer_dataframe(
    clean_seq: List[SeqRecord],
    melting_temperature: float,
    polymerase: str,
    primer_concentration: float,
    up_homology: str,
    dw_homology: str,
    min_primer_len=10,
) -> pd.DataFrame:
    """
    Generates a DataFrame containing primer sequences, melting temperatures,
    annealing temperature, and incremental primer names, appending specified homology
    sequences to each primer.

    Parameters
    ----------
    clean_seq : list of SeqRecord
        List of clean sequences for which to generate primers.
    melting_temperature : float
        Target melting temperature for the primers.
    polymerase : str
        Polymerase used for primer Tm calculations.
    primer_concentration : float
        Concentration of primers used in the reaction.
    up_homology : str
        Upstream homology sequence to append to forward primers.
    dw_homology : str
        Downstream homology sequence to append to reverse primers.
    min_primer_len : int
        Minimum length of primers - starts at this lenght and then finds optimal mt.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the columns: template, f_primer_anneal(5-3),
        r_primer_anneal(5-3), f_tm, r_tm, ta, f_primer_sequences(5-3),
        r_primer_sequences(5-3), f_primer_name, r_primer_name.

    Examples
    --------
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> clean_seq = [SeqRecord(Seq("ATGC"), id="gene1"), SeqRecord(Seq("GCAT"), id="gene2")]
    >>> df = generate_primer_dataframe(clean_seq, 60, 'Taq', 0.5, 'ATCG', 'CGTA', 1)
    >>> print(df.head())
    """
    amplicons = make_amplicons(
        clean_seq,
        target_tm=melting_temperature,
        limit=min_primer_len,
        polymerase=polymerase,
        primer_concentration=primer_concentration,
    )
    data = []

    for amplicon, gene_id in zip(amplicons, clean_seq):
        forward_primer_seq = str(amplicon.forward_primer.seq)
        reverse_primer_seq = str(amplicon.reverse_primer.seq)
        f_tm = primer_tm_neb(
            forward_primer_seq, prodcode=polymerase, conc=primer_concentration
        )
        r_tm = primer_tm_neb(
            reverse_primer_seq, prodcode=polymerase, conc=primer_concentration
        )
        ta = primer_ta_neb(
            forward_primer_seq,
            reverse_primer_seq,
            prodcode=polymerase,
            conc=primer_concentration,
        )

        data.append(
            {
                "template": f"{str(gene_id.id)}",
                "f_primer_anneal(5-3)": forward_primer_seq,
                "r_primer_anneal(5-3)": reverse_primer_seq,
                "f_tm": f_tm,
                "r_tm": r_tm,
                "ta": ta,
            }
        )

    df = pd.DataFrame(data)

    df["f_primer_sequences(5-3)"] = up_homology + df["f_primer_anneal(5-3)"]
    df["r_primer_sequences(5-3)"] = dw_homology + df["r_primer_anneal(5-3)"]

    f_primer_names = [f"primer_fwd_{df.loc[i, 'template']}" for i in range(len(df))]
    r_primer_names = [f"primer_rev_{df.loc[i, 'template']}" for i in range(len(df))]

    df["f_primer_name"] = f_primer_names
    df["r_primer_name"] = r_primer_names

    return df


# TODO merge with the function below and reorder primers
def create_idt_order_dataframe(
    df: pd.DataFrame, concentration: str = "25nm", purification: str = "STD"
) -> pd.DataFrame:
    """
    Creates a DataFrame formatted for IDT (Integrated DNA Technologies) orders,
    containing primer names, sequences, specified concentration, and purification method.
    Allows for customization of the concentration and purification while providing
    default values.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the primer information. Expected to have columns
        'f_primer_name', 'r_primer_name', 'f_primer_sequences(5-3)', and 'r_primer_sequences(5-3)'.
    concentration : str, optional
        The concentration for all primers, by default "25nm".
    purification : str, optional
        The purification method for all primers, by default "STD".

    Returns
    -------
    pd.DataFrame
        A new DataFrame with columns 'Name', 'Sequences', 'Concentration', and 'Purification',
        formatted for IDT orders, with the specified concentration and purification methods.

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     'f_primer_name': ['f_primer1', 'f_primer2'],
    ...     'r_primer_name': ['r_primer1', 'r_primer2'],
    ...     'f_primer_sequences(5-3)': ['ATCG', 'CGTA'],
    ...     'r_primer_sequences(5-3)': ['GATC', 'TACG']
    ... })
    >>> idt_df = create_idt_order_dataframe(df, concentration="100nm", purification="PAGE")
    >>> print(idt_df)
    """
    primer_names = df["f_primer_name"].tolist() + df["r_primer_name"].tolist()
    primer_sequences = (
        df["f_primer_sequences(5-3)"].tolist() + df["r_primer_sequences(5-3)"].tolist()
    )

    idt_data = {
        "Name": primer_names,
        "Sequence": primer_sequences,
        "Concentration": [concentration] * len(primer_names),
        "Purification": [purification] * len(primer_names),
    }

    idt_df = pd.DataFrame(idt_data)

    return idt_df


# TODO almost exactly like the funciton above
def primers_to_IDT(
    sgRNAs_p: list, concentration: str = "25nm", purification: str = "STD"
) -> pd.DataFrame:
    """ "Generates a DataFrame containing primers that
    are ready for submission to IDT for synthesis,

    Parameters
    ----------
    sgRNAs_p : list
        list of pydna.Dseqrecords
    concentration : str
    purification : str

    Returns
    --------
    best_gRNAs : pd.DataFrame
        Top sgRNAs

    """

    name_lst = [name.name for name in sgRNAs_p]
    seq_lst = [str(seq.seq) for seq in sgRNAs_p]

    d = {
        "Name": name_lst,
        "Sequence": seq_lst,
        "Concentration": concentration,
        "Purification": purification,
    }
    df = pd.DataFrame(d)
    return df


# TODO rename and merge with the dataframe_to_seqrecords funciton in utils.py
def make_primer_records(filtered_df: pd.DataFrame) -> List[Dseqrecord]:
    """
    Create Dseqrecord objects for the forward and reverse primers with names according to gene and location.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        DataFrame containing the filtered spacer sequences with columns for forward and reverse primers.

    Returns
    -------
    List[Dseqrecord]
        A list of Dseqrecord objects for the forward and reverse primers.
    """
    primer_records: List[Dseqrecord] = []

    for index, row in filtered_df.iterrows():
        # Create Dseqrecord for forward primer
        fwd_primer_record = Dseqrecord(
            row["Fwd Primer"],
            id=f"{row['locus_tag']}_fwd_{str(index)}",
            name=f"{row['locus_tag']}_loc_{row['sgrna_loc']}_fwd",
        )
        fwd_label = f"Fwd_{row['locus_tag']}_loc_{row['sgrna_loc']}"
        fwd_primer_record.features.append(
            SeqFeature(
                FeatureLocation(0, len(row["Fwd Primer"])),
                type="primer_bind",
                qualifiers={"label": fwd_label},
            )
        )

        # Create Dseqrecord for reverse primer
        rev_primer_record = Dseqrecord(
            row["Rev Primer"],
            id=f"{row['locus_tag']}_rev_{str(index)}",
            name=f"{row['locus_tag']}_loc_{row['sgrna_loc']}_rev",
        )
        rev_label = f"Rev_{row['locus_tag']}_loc_{row['sgrna_loc']}"
        rev_primer_record.features.append(
            SeqFeature(
                FeatureLocation(0, len(row["Rev Primer"])),
                type="primer_bind",
                qualifiers={"label": rev_label},
            )
        )

        primer_records.extend([fwd_primer_record, rev_primer_record])

    return primer_records


def checking_primers(
    record,
    locus_tags,
    flanking_region=500,
    target_tm=65,
    limit=10,
    primer_concentration=0.4,
    polymerase="onetaq-3",
    **primer_kwargs,
):
    if not primer_kwargs:
        primer_kwargs = {"conc": primer_concentration, "prodcode": polymerase}

    primer_info = []

    for locus_tag in locus_tags:
        feature = None
        for f in record.features:
            if (
                f.type == "CDS"
                and "locus_tag" in f.qualifiers
                and f.qualifiers["locus_tag"][0] == locus_tag
            ):
                feature = f
                break

        if feature is None:
            raise ValueError(f"Locus tag {locus_tag} not found in the GenBank file.")

        start = max(0, feature.location.start - flanking_region)
        end = min(len(record), feature.location.end + flanking_region)
        target_seq = record.seq[start:end]

        dseqrecord = Dseqrecord(target_seq)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            primers = primer_design(
                dseqrecord, target_tm=target_tm, limit=limit, **primer_kwargs
            )

            # Check for the specific warning about non-unique PCR products
            for warning in w:
                if "unique PCR product" in str(warning.message):
                    raise ValueError(
                        "Primers do not yield a unique PCR product, iterating again."
                    )

        forward_primer_name = f"{locus_tag}_fwd_checking_primer"
        reverse_primer_name = f"{locus_tag}_rev_checking_primer"
        primers.forward_primer.name = forward_primer_name
        primers.reverse_primer.name = reverse_primer_name

        primer_info.append(
            {
                "locus tag": locus_tag,
                "f_primer_name": forward_primer_name,
                "r_primer_name": reverse_primer_name,
                "f_primer_sequences(5-3)": str(primers.forward_primer.seq),
                "r_primer_sequences(5-3)": str(primers.reverse_primer.seq),
                "f_tm": primer_tm_neb(str(primers.forward_primer.seq), **primer_kwargs),
                "r_tm": primer_tm_neb(str(primers.reverse_primer.seq), **primer_kwargs),
                "ta": primer_ta_neb(
                    str(primers.forward_primer.seq),
                    str(primers.reverse_primer.seq),
                    **primer_kwargs,
                ),
            }
        )

    primer_df = pd.DataFrame(primer_info)

    return primer_df


def validate_primers(analysis_df):
    for _, row in analysis_df.iterrows():
        if (
            row["homodimer_forward_deltaG (kcal/mol)"] < -9
            or row["homodimer_reverse_deltaG (kcal/mol)"] < -9
            or row["heterodimer_deltaG (kcal/mol)"] < -9
            or row["hairpin_forward_deltaG (kcal/mol)"] < -9
            or row["hairpin_reverse_deltaG (kcal/mol)"] < -9
            or row["hairpin_forward_structure_found"]
            or row["hairpin_reverse_structure_found"]
            or row["f_tm"] >= 72
            or row["r_tm"] >= 72
            or abs(row["f_tm"] - row["r_tm"]) > 4
        ):
            return False
    return True


def find_best_check_primers_from_genome(
    record,
    locus_tags,
    flanking_region=500,
    target_tm=65,
    limit=10,
    primer_concentration=0.4,
    polymerase="onetaq-3",
    max_iterations=50,
    **primer_kwargs,
):
    if not primer_kwargs:
        primer_kwargs = {"conc": primer_concentration, "prodcode": polymerase}

    initial_flanking_region = flanking_region
    successful_primers = []
    remaining_locus_tags = set(locus_tags.copy())

    for iteration in range(max_iterations):
        if not remaining_locus_tags:
            break

        try:
            primer_df = checking_primers(
                record=record,
                locus_tags=list(remaining_locus_tags),
                flanking_region=flanking_region,
                target_tm=target_tm,
                limit=limit,
                primer_concentration=primer_concentration,
                polymerase=polymerase,
                **primer_kwargs,
            )

            analysis_df = analyze_primers_and_hairpins(primer_df)

            valid_primers = []
            invalid_primers = set()

            for idx, row in primer_df.iterrows():
                locus_tag = row["locus tag"]
                if validate_primers(analysis_df.iloc[[idx]]):
                    row["flanking_region"] = flanking_region
                    valid_primers.append(
                        pd.concat(
                            [
                                row.to_frame().T.reset_index(drop=True),
                                analysis_df.iloc[[idx]].reset_index(drop=True),
                            ],
                            axis=1,
                        )
                    )
                    remaining_locus_tags.discard(
                        locus_tag
                    )  # Remove successful locus tag
                else:
                    invalid_primers.add(locus_tag)

            successful_primers.extend(valid_primers)
            remaining_locus_tags = invalid_primers

            if not valid_primers:
                flanking_region += (
                    1  # Increase flanking region only if no valid primers were found
                )

        except ValueError as e:
            if "unique PCR product" in str(e):
                flanking_region += 50
            else:
                raise

    if not successful_primers:
        raise ValueError(
            f"Could not find suitable primers within the maximum number of iterations, started with {initial_flanking_region} flanking region and ended with {flanking_region} flanking region"
        )

    final_result_df = pd.concat(successful_primers, axis=0).reset_index(drop=True)
    return final_result_df
