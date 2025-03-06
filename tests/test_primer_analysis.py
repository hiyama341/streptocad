import pandas as pd
import pytest

from streptocad.primers.primer_analysis import analyze_primers_and_hairpins


@pytest.fixture
def primer_df():
    primer_df = pd.read_csv("tests/test_files/primer_df_workflow1.csv")
    return primer_df


@pytest.fixture
def expected_primer_analysis():
    primer_df = pd.read_csv("tests/test_files/expected_primer_analysis_workflow1.csv")
    return primer_df


def test_analyze_primers_and_hairpins_basic(primer_df, expected_primer_analysis):
    """
    Test the analyze_primers_and_hairpins function with basic input
    to ensure that it produces the correct output.
    """
    result_df = analyze_primers_and_hairpins(primer_df)

    # Print columns to identify discrepancies
    print("Result DataFrame columns:", result_df.columns)
    print("Expected DataFrame columns:", expected_primer_analysis.columns)

    # Check for column mismatches
    result_cols = set(result_df.columns)
    expected_cols = set(expected_primer_analysis.columns)
    if result_cols != expected_cols:
        print("Columns missing in result_df:", expected_cols - result_cols)
        print(
            "Columns missing in expected_primer_analysis:", result_cols - expected_cols
        )

    # Align both DataFrames to have the same columns
    common_cols = list(result_cols.intersection(expected_cols))
    result_df = result_df[common_cols]
    expected_primer_analysis = expected_primer_analysis[common_cols]

    # Rounding numerical values for comparison due to potential floating-point discrepancies
    numeric_cols = [
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
    result_df.loc[:, numeric_cols] = result_df[numeric_cols].round(2)
    expected_primer_analysis.loc[:, numeric_cols] = expected_primer_analysis[
        numeric_cols
    ].round(2)

    # Assert that the result DataFrame matches the expected DataFrame
    pd.testing.assert_frame_equal(result_df, expected_primer_analysis, check_like=True)


def test_analyze_primers_and_hairpins_edge_cases():
    """
    Test that the function correctly skips primers longer than 60 nucleotides.
    """
    edge_case_df = pd.DataFrame(
        {
            "f_primer_sequences(5-3)": ["A" * 61, "TGCATGCATGC"],
            "r_primer_sequences(5-3)": ["G" * 61, "CGTACGTACGT"],
            "f_primer_name": ["forward_long", "forward_valid"],
            "r_primer_name": ["reverse_long", "reverse_valid"],
            "ta": [60.0, 60.0],
        }
    )

    result_df = analyze_primers_and_hairpins(edge_case_df)

    # Ensure only the valid pair is processed
    assert len(result_df) == 1
    assert result_df.iloc[0]["primer_pair"] == "forward_valid & reverse_valid"
