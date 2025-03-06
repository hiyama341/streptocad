import pandas as pd
import pytest
from pydna.dseqrecord import Dseqrecord
from Bio import SeqIO


from streptocad.primers.primer_generation import (
    generate_primer_dataframe,
    create_idt_order_dataframe,
    primers_to_IDT,
    make_primer_records,
    checking_primers,
    validate_primers,
    find_best_check_primers_from_genome,
)
from streptocad.sequence_loading.sequence_loading import (
    load_and_process_genome_sequences,
)
from streptocad.utils import polymerase_dict


@pytest.fixture
def primer_df():
    primer_df = pd.read_csv("tests/test_files/primer_df_workflow1.csv")
    return primer_df


@pytest.fixture
def expected_primer_analysis():
    primer_df = pd.read_csv("tests/test_files/expected_primer_analysis_workflow1.csv")
    return primer_df


@pytest.fixture
def clean_sequences():
    path_to_file = "tests/test_files/GOE_regulators.gb"
    clean_seq = load_and_process_genome_sequences(path_to_file)

    return clean_seq


@pytest.fixture
def sgRNAs_p():
    # Creating a mock list of Dseqrecord objects as a fixture
    primer1 = Dseqrecord("TCCACCGGCGCCGCGTCCAG", name="primer_1")
    primer2 = Dseqrecord("GCTGGGCGCGATCGGCTCGC", name="primer_2")
    return [primer1, primer2]


@pytest.fixture
def coelicolor_genbank_record():
    filepath = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"
    record = Dseqrecord(SeqIO.read(filepath, "genbank"))
    return record


@pytest.fixture
def filtered_df():
    data = {
        "strain_name": ["NC_003888.3", "NC_003888.3", "NC_003888.3", "NC_003888.3"],
        "locus_tag": ["SCO5087", "SCO5087", "SCO5090", "SCO5090"],
        "gene_loc": [5529801, 5529801, 5532706, 5532706],
        "gene_strand": [1, 1, 1, 1],
        "sgrna_strand": [1, -1, -1, -1],
        "sgrna_loc": [131, 163, 199, 505],
        "gc": [0.558824, 0.617647, 0.764706, 0.735294],
        "pam": ["TTC", "TTC", "TTC", "TTC"],
        "sgrna": [
            "AATCCGAACAGCTCCTTCGATGCGTTCGTCCGGT",
            "GGATTGAAGCGCAGAGTCGTCATCACGGGCGTCG",
            "CGGATCTGGGCGCGCGTGGGCGGCCGGGTGAAGA",
            "ACGTTCGAGGACACGCTCCGGGTGCCGTCCGGGG",
        ],
        "sgrna_seed_sequence": ["AATCCGAA", "GGATTGAA", "CGGATCTG", "ACGTTCGA"],
        "off_target_count": [0, 0, 3, 4],
        "Fwd Primer": [
            "AATCCGAACAGCTCCTTCGATGCGTTCGTCCGGTGTCGCcCggCaaAaccGg",
            "GGATTGAAGCGCAGAGTCGTCATCACGGGCGTCGGTCGCcCggCaaAaccGg",
            "CGGATCTGGGCGCGCGTGGGCGGCCGGGTGAAGAGTCGCcCggCaaAaccGg",
            "ACGTTCGAGGACACGCTCCGGGTGCCGTCCGGGGGTCGCcCggCaaAaccGg",
        ],
        "Rev Primer": [
            "ACCGGACGAACGCATCGAAGGAGCTGTTCGGATTGTTTCAATCCACGCGCCCGT",
            "CGACGCCCGTGATGACGACTCTGCGCTTCAATCCGTTTCAATCCACGCGCCCGT",
            "TCTTCACCCGGCCGCCCACGCGCGCCCAGATCCGGTTTCAATCCACGCGCCCGT",
            "CCCCGGACGGCACCCGGAGCGTGTCCTCGAACGTGTTTCAATCCACGCGCCCGT",
        ],
    }

    filtered_df = pd.DataFrame(data)

    return filtered_df


@pytest.fixture
def find_best_checking_primers_df():
    data = {
        "locus tag": ["SCO5087"],
        "f_primer_name": ["SCO5087_fwd_checking_primer"],
        "r_primer_name": ["SCO5087_rev_checking_primer"],
        "f_primer_sequences(5-3)": ["GACGATTCGGCCCGTG"],
        "r_primer_sequences(5-3)": ["CAGGGCGTCCAGGC"],
        "f_tm": [59],
        "r_tm": [58],
        "ta": [61],
        "flanking_region": [500],
        "annealing_temperature": [61],
        "primer_pair": ["SCO5087_fwd_checking_primer & SCO5087_rev_checking_primer"],
        "homodimer_forward_tm": [9.596675850895167],
        "homodimer_forward_deltaG (kcal/mol)": [-0.3393059724374107],
        "homodimer_reverse_tm": [9.70624571345519],
        "homodimer_reverse_deltaG (kcal/mol)": [-3.4154272069089346],
        "heterodimer_tm": [1.9185669591946635],
        "heterodimer_deltaG (kcal/mol)": [-2.396641324148026],
        "hairpin_forward_structure_found": [False],
        "hairpin_forward_tm": [0.0],
        "hairpin_forward_deltaG (kcal/mol)": [0.0],
        "hairpin_reverse_structure_found": [False],
        "hairpin_reverse_tm": [0.0],
        "hairpin_reverse_deltaG (kcal/mol)": [0.0],
    }

    df = pd.DataFrame(data)

    return df


@pytest.fixture
def checking_primers_df():
    data = {
        "locus tag": ["SCO5087"],
        "f_primer_name": ["SCO5087_fwd_checking_primer"],
        "r_primer_name": ["SCO5087_rev_checking_primer"],
        "f_primer_sequences(5-3)": ["GACGATTCGGCCCGTG"],
        "r_primer_sequences(5-3)": ["CAGGGCGTCCAGGC"],
        "f_tm": [59],
        "r_tm": [58],
        "ta": [61],
    }

    df = pd.DataFrame(data)

    return df


def test_generate_primer_dataframe_with_fixture(primer_df, clean_sequences):
    # Define the parameters for the test
    # 3 Choose overlapping sequences for our plasmid we can use the following
    up_homology = "GGCGAGCAACGGAGGTACGGACAGG".upper()
    dw_homology = "CGCAAGCCGCCACTCGAACGGAAGG".upper()

    #### Advanced settings ####
    # 4 Choose polymerase and target melting temperature
    chosen_polymerase = "Phusion High-Fidelity DNA Polymerase (GC Buffer)"
    # chosen_polymerase = 'default'

    melting_temperature = 60
    primer_concentration = 0.4
    # Call the function using the clean_sequences fixture
    result_df = generate_primer_dataframe(
        clean_seq=clean_sequences,
        melting_temperature=melting_temperature,
        polymerase=polymerase_dict[chosen_polymerase],
        primer_concentration=primer_concentration,
        up_homology=up_homology,
        dw_homology=dw_homology,
        min_primer_len=15,
    )

    # Prepare primer_df to match expected columns for comparison
    expected_df = primer_df.copy()
    expected_df["f_tm"] = expected_df["f_tm"].round(2)
    expected_df["r_tm"] = expected_df["r_tm"].round(2)
    expected_df["ta"] = expected_df["ta"].round(2)

    # Check only relevant columns for comparison (ignore other potential columns in primer_df)
    expected_columns = [
        "template",
        "f_primer_anneal(5-3)",
        "r_primer_anneal(5-3)",
        "f_tm",
        "r_tm",
        "ta",
        "f_primer_sequences(5-3)",
        "r_primer_sequences(5-3)",
        "f_primer_name",
        "r_primer_name",
    ]
    result_df = result_df[expected_columns]
    expected_df = expected_df[expected_columns]

    # Round numerical values for comparison to avoid floating-point discrepancies
    numeric_cols = ["f_tm", "r_tm", "ta"]
    result_df[numeric_cols] = result_df[numeric_cols].round(2)
    expected_df[numeric_cols] = expected_df[numeric_cols].round(2)

    # Assert that the result DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(result_df, expected_df, check_like=True)


def test_create_idt_order_dataframe(primer_df):
    # Define the parameters for the test
    concentration = "100nm"
    purification = "PAGE"

    # Call the function using the primer_df fixture
    idt_df = create_idt_order_dataframe(
        primer_df, concentration=concentration, purification=purification
    )

    # Prepare the expected DataFrame
    expected_names = (
        primer_df["f_primer_name"].tolist() + primer_df["r_primer_name"].tolist()
    )
    expected_sequences = (
        primer_df["f_primer_sequences(5-3)"].tolist()
        + primer_df["r_primer_sequences(5-3)"].tolist()
    )
    expected_concentration = [concentration] * len(expected_names)
    expected_purification = [purification] * len(expected_names)

    expected_data = {
        "Name": expected_names,
        "Sequence": expected_sequences,
        "Concentration": expected_concentration,
        "Purification": expected_purification,
    }
    expected_df = pd.DataFrame(expected_data)

    # Assert that the result DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(idt_df, expected_df)


def test_primers_to_IDT(sgRNAs_p):
    # Define parameters for the test
    concentration = "100nm"
    purification = "PAGE"

    # Call the function using the sgRNAs_p fixture
    idt_df = primers_to_IDT(
        sgRNAs_p, concentration=concentration, purification=purification
    )

    # Prepare the expected DataFrame
    expected_data = {
        "Name": ["primer_1", "primer_2"],
        "Sequence": ["TCCACCGGCGCCGCGTCCAG", "GCTGGGCGCGATCGGCTCGC"],
        "Concentration": [concentration] * 2,
        "Purification": [purification] * 2,
    }
    expected_df = pd.DataFrame(expected_data)

    # Assert that the result DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(idt_df, expected_df)


def test_make_primer_records(filtered_df):
    # Call the function with the test fixture
    primer_records = make_primer_records(filtered_df)

    # Check that the length of the returned list is twice the number of rows in the input DataFrame
    assert len(primer_records) == len(filtered_df) * 2, (
        "The number of primer records should be twice the number of DataFrame rows."
    )

    # Check each record to see if it matches the expected output
    for i, record in enumerate(primer_records):
        row = filtered_df.iloc[i // 2]
        primer_type = "fwd" if i % 2 == 0 else "rev"
        primer_seq = row["Fwd Primer"] if primer_type == "fwd" else row["Rev Primer"]

        # Verify sequence
        assert str(record.seq) == primer_seq, (
            f"Primer sequence does not match for record {record.id}."
        )

        # Verify ID and Name
        expected_id = f"{row['locus_tag']}_{primer_type}_{i // 2}"
        expected_name = f"{row['locus_tag']}_loc_{row['sgrna_loc']}_{primer_type}"
        assert record.id == expected_id, (
            f"Primer ID does not match for record {record.id}."
        )
        assert record.name == expected_name, (
            f"Primer name does not match for record {record.id}."
        )

        # Calculate expected feature label
        expected_feature_label = (
            f"{primer_type.capitalize()}_{row['locus_tag']}_loc_{row['sgrna_loc']}"
        )

        # Debug: Print feature labels
        feature_labels = [
            feat.qualifiers["label"][0]
            for feat in record.features
            if "label" in feat.qualifiers
        ]
        print(
            f"Expected feature label: {expected_feature_label}, Actual feature labels: {feature_labels}"
        )


def test_checking_primers(coelicolor_genbank_record, checking_primers_df):
    # Define the locus tags to check
    locus_tags = ["SCO5087"]

    # Call the checking_primers function
    result_df = checking_primers(
        coelicolor_genbank_record,
        locus_tags,
        flanking_region=500,
        target_tm=65,
        limit=10,
        primer_concentration=0.4,
        polymerase="phusion-1",
    )

    # Define expected numeric columns
    numeric_cols = ["f_tm", "r_tm", "ta"]

    # Round numeric columns for comparison
    result_df[numeric_cols] = result_df[numeric_cols].round(2)
    checking_primers_df[numeric_cols] = checking_primers_df[numeric_cols].round(2)

    # Check for missing columns in result_df and add them if necessary
    for col in checking_primers_df.columns:
        if col not in result_df.columns:
            result_df[col] = None

    # Reorder columns to match the expected DataFrame
    result_df = result_df[checking_primers_df.columns]

    # Assert that the resulting DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(result_df, checking_primers_df, check_like=True)


# Sample test data for validate_primers function
@pytest.fixture
def analysis_df_valid():
    # Data simulates a case where all primers pass the validation
    data = {
        "homodimer_forward_deltaG (kcal/mol)": [-8.5, -8.0],
        "homodimer_reverse_deltaG (kcal/mol)": [-8.0, -8.5],
        "heterodimer_deltaG (kcal/mol)": [-7.5, -8.0],
        "hairpin_forward_deltaG (kcal/mol)": [-5.0, -4.5],
        "hairpin_reverse_deltaG (kcal/mol)": [-4.5, -5.0],
        "hairpin_forward_structure_found": [False, False],
        "hairpin_reverse_structure_found": [False, False],
    }
    return pd.DataFrame(data)


@pytest.fixture
def analysis_df_invalid():
    # Data simulates a case where some primers fail the validation
    data = {
        "homodimer_forward_deltaG (kcal/mol)": [-10.0, -8.0],
        "homodimer_reverse_deltaG (kcal/mol)": [-8.0, -10.0],
        "heterodimer_deltaG (kcal/mol)": [-7.5, -8.0],
        "hairpin_forward_deltaG (kcal/mol)": [-5.0, -10.5],
        "hairpin_reverse_deltaG (kcal/mol)": [-4.5, -5.0],
        "hairpin_forward_structure_found": [False, True],
        "hairpin_reverse_structure_found": [False, False],
    }
    return pd.DataFrame(data)


def test_validate_primers_valid(analysis_df_valid):
    # Validate primers should return True as all conditions are met
    assert validate_primers(analysis_df_valid) == True


def test_validate_primers_invalid(analysis_df_invalid):
    # Validate primers should return False as some conditions are not met
    assert validate_primers(analysis_df_invalid) == False


def test_find_best_check_primers_from_genome(
    coelicolor_genbank_record, find_best_checking_primers_df
):
    # Define the locus tags to check
    locus_tags = ["SCO5087"]

    # Call the find_best_check_primers_from_genome function
    result_df = find_best_check_primers_from_genome(
        record=coelicolor_genbank_record,
        locus_tags=locus_tags,
        flanking_region=500,
        target_tm=65,
        limit=10,
        primer_concentration=0.4,
        polymerase="phusion-1",
        max_iterations=50,
    )

    # Ensure both DataFrames have the same data types for comparison
    numeric_cols = [
        "f_tm",
        "r_tm",
        "ta",
        "homodimer_forward_tm",
        "homodimer_forward_deltaG (kcal/mol)",
        "homodimer_reverse_tm",
        "homodimer_reverse_deltaG (kcal/mol)",
        "heterodimer_tm",
        "heterodimer_deltaG (kcal/mol)",
        "hairpin_forward_tm",
        "hairpin_forward_deltaG (kcal/mol)",
        "hairpin_reverse_tm",
        "hairpin_reverse_deltaG (kcal/mol)",
    ]

    # Convert numeric columns to float for consistency
    result_df[numeric_cols] = result_df[numeric_cols].astype(float)
    find_best_checking_primers_df[numeric_cols] = find_best_checking_primers_df[
        numeric_cols
    ].astype(float)

    # Convert the flanking_region column to int for consistency
    result_df["flanking_region"] = result_df["flanking_region"].astype(int)
    find_best_checking_primers_df["flanking_region"] = find_best_checking_primers_df[
        "flanking_region"
    ].astype(int)

    # Round the floating point numbers to match the expected precision
    result_df[numeric_cols] = result_df[numeric_cols].round(6)
    find_best_checking_primers_df[numeric_cols] = find_best_checking_primers_df[
        numeric_cols
    ].round(6)

    # Assert that the result DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(
        result_df, find_best_checking_primers_df, check_like=True
    )
