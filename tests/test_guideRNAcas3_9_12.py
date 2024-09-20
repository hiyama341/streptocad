from streptocad.crispr.guideRNAcas3_9_12 import (
    SgRNAargs, 
    revcomp, 
    parse_genbank_record,
    find_off_target_hits, 
    find_sgrna_hits_cas9,
    find_sgrna_hits_cas12a, 
    find_sgrna_hits_cas3, 
    filter_guides,
    extract_sgRNAs,
)

import pytest
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd

import pytest
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

@pytest.fixture
def coelicolor_genbank_record():
    filepath = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"
    record = Dseqrecord(SeqIO.read(filepath, "genbank"))
    return record



### TESTS ###
def test_sgrnaargs_initialization():
    # Create a mock Dseqrecord object
    seq = Seq("ATGCATGCATGCATGCATGC")
    dseqrecord = Dseqrecord(seq, linear=True)
    dseqrecord.id = "Mock_Strain"

    # Valid initialization
    args = SgRNAargs(
        dseqrecord=dseqrecord,
        locus_tag=['gene1', 'gene2'],
        cas_type='cas9',
        downstream=20,
        off_target_seed=10,
        gc_upper=0.8,
        gc_lower=0.2,
        pam_remove=['PAM1'],
        downstream_remove=['remove1'],
        sgrna_remove=['sgrna1'],
        off_target_upper=12,
        extension_to_promoter_region=5,
        upstream_tss=150,
        dwstream_tss=150,
        target_non_template_strand=False
    )

    # Test attributes
    assert args.dseqrecord == dseqrecord
    assert args.locus_tag == ['gene1', 'gene2']
    assert args.cas_type == 'cas9'
    assert args.downstream == 20
    assert args.off_target_seed == 10
    assert args.gc_upper == 0.8
    assert args.gc_lower == 0.2
    assert args.pam_remove == ['PAM1']
    assert args.downstream_remove == ['remove1']
    assert args.sgrna_remove == ['sgrna1']
    assert args.off_target_upper == 12
    assert args.extension_to_promoter_region == 5
    assert args.upstream_tss == 150
    assert args.dwstream_tss == 150
    assert args.target_non_template_strand == False
    assert args.strain_name == "Mock_Strain"

def test_invalid_dseqrecord_initialization():
    # Invalid initialization: dseqrecord is not a Dseqrecord instance
    with pytest.raises(ValueError):
        SgRNAargs(
            dseqrecord="Not a Dseqrecord",  # Invalid type
            locus_tag=['gene1'],
        )

def test_default_values():
    # Create a mock Dseqrecord object
    seq = Seq("ATGCATGCATGCATGCATGC")
    dseqrecord = Dseqrecord(seq, linear=True)
    dseqrecord.id = "Mock_Strain"

    # Test initialization with defaults
    args = SgRNAargs(
        dseqrecord=dseqrecord,
        locus_tag=['gene1']
    )

    # Check default values
    assert args.step == ['find', 'filter']
    assert args.downstream == 15
    assert args.off_target_seed == 9
    assert args.pam_remove is None
    assert args.downstream_remove is None
    assert args.sgrna_remove is None
    assert args.gc_upper == 1
    assert args.gc_lower == 0
    assert args.off_target_upper == 10
    assert args.extension_to_promoter_region == 0
    assert args.upstream_tss == 100
    assert args.dwstream_tss == 100
    assert args.target_non_template_strand == True


def test_revcomp_basic():
    # Basic sequences
    assert revcomp("ATGC") == "GCAT"
    assert revcomp("AAGCTT") == "AAGCTT"[::-1].translate(str.maketrans("ATGC", "TACG"))  # "AAGCTT" -> "AAGCTT" -> "TTCGAA"
    assert revcomp("CGTA") == "TACG"

def test_revcomp_edge_cases():
    # Edge cases
    assert revcomp("") == ""  # Empty string
    assert revcomp("A") == "T"  # Single nucleotide
    assert revcomp("C") == "G"
    assert revcomp("G") == "C"
    assert revcomp("T") == "A"


def test_parse_genbank_record_basic(coelicolor_genbank_record):
    # Call the function with the actual GenBank record
    sequences = parse_genbank_record(coelicolor_genbank_record)

    # Debug print to check sequences and DataFrame
    print("Sequences:", sequences)

    # Test sequences (example check, update according to what you expect)
    assert len(sequences) == 2  # Check the length of the seq


def test_parse_genbank_record_no_genes(coelicolor_genbank_record):
    # Modify the record to have no gene features (if applicable)
    coelicolor_genbank_record.features = []

    # Call the function
    sequences = parse_genbank_record(coelicolor_genbank_record)

    # Test sequences
    assert len(sequences) == 2  # Check that both strands are returned
    assert len(sequences[0]) == 8667507  # Ensure the watson sequence 
    assert len(sequences[1]) == 8667507  # Ensure the reverse complement is not empty


def test_parse_genbank_record_partial_info(coelicolor_genbank_record):
    # Modify the record to have partial gene information
    for feature in coelicolor_genbank_record.features:
        if feature.type == "gene":
            feature.qualifiers.pop("gene", None)  # Remove the gene name if present

    # Call the function
    sequences = parse_genbank_record(coelicolor_genbank_record)

    # Test sequences
    assert len(sequences) == 2  # Check that both strands are returned
    assert len(sequences[0]) == 8667507  # Ensure the watson sequence 
    assert len(sequences[1]) == 8667507  # Ensure the reverse complement is not empty



@pytest.fixture
def sgrna_args_cas9(coelicolor_genbank_record):
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=['SCO5087'],  
        cas_type='cas9',
        gc_upper=0.999,
        gc_lower=0.01,
        off_target_seed=13,
        off_target_upper=10,
        step=['find', 'filter'] 
    )


@pytest.fixture
def sgrna_args_cas3(coelicolor_genbank_record):
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=['SCO5087','SCO5089', 'SCO5090'],  
        cas_type='cas3',
        gc_upper=0.8,
        gc_lower=0.2,
        off_target_seed=8,
        off_target_upper=10,
        step=['find', 'filter'] 
    )

@pytest.fixture
def expected_sgrna_df():
    filepath = "tests/test_files/unfiltered_sgrna_SCO5087_df.csv"
    df = pd.read_csv(filepath, index_col=0)  # Load the file, treating the first column as the index
    df.reset_index(drop=True, inplace=True)  # Reset the index to remove the index column
    return df

@pytest.fixture
def expected_cas3_sgrna_df():
    filepath = "tests/test_files/unfiltered_sgrna_SCO5087_SCO5089_SCO5090_cas3_df.csv"
    df = pd.read_csv(filepath, index_col=0) 
    df.reset_index(drop=True, inplace=True) 
    return df


def test_find_off_target_hits_cas9(coelicolor_genbank_record, sgrna_args_cas9):
    # Extract sequences from the Dseqrecord object
    sequences = [str(coelicolor_genbank_record.seq), str(coelicolor_genbank_record.reverse_complement().seq)]

    # Use the function to find off-target hits
    off_target_counter = find_off_target_hits(sequences, sgrna_args_cas9.off_target_seed, cas_type='cas9')
    print(off_target_counter)
    # Check that the Counter is not empty
    assert off_target_counter, "The off-target Counter should not be empty"

    # Check if certain expected seed sequences are present
    # (Here, you would check against known or expected off-target sequences for the Cas9 system)
    expected_seeds = ['GCGCCGCGTCCAG', 'GCGATCGGCTCGC']  # Replace with actual expected sequences if known
    for seed in expected_seeds:
        assert seed in off_target_counter, f"Expected seed sequence '{seed}' not found in off-target hits"

    # Check that the counts are integers and make sense
    for count in off_target_counter.values():
        assert isinstance(count, int), "The count for off-targets should be an integer"



def test_find_off_target_hits_cas12a(coelicolor_genbank_record):
    # Extract sequences from the Dseqrecord object
    sequences = [str(coelicolor_genbank_record.seq), str(coelicolor_genbank_record.reverse_complement().seq)]

    # Use the function to find off-target hits for Cas12a
    off_target_counter = find_off_target_hits(sequences, off_target_seed=7, cas_type='cas12a')

    # Check that the Counter is not empty
    assert off_target_counter, "The off-target Counter should not be empty for Cas12a"

    # Check if certain expected seed sequences are present for Cas12a
    expected_seeds = ['ACGTGAA', 'TTCGGCA']  # Replace with actual expected sequences if known
    for seed in expected_seeds:
        assert seed in off_target_counter, f"Expected seed sequence '{seed}' not found in off-target hits for Cas12a"



def test_find_off_target_hits_cas3(coelicolor_genbank_record):
    # Extract sequences from the Dseqrecord object
    sequences = [str(coelicolor_genbank_record.seq), str(coelicolor_genbank_record.reverse_complement().seq)]

    # Use the function to find off-target hits for Cas3
    off_target_counter = find_off_target_hits(sequences, off_target_seed=10, cas_type='cas3')

    # Check that the Counter is not empty
    assert off_target_counter, "The off-target Counter should not be empty for Cas3"

    # Check if certain expected seed sequences are present for Cas3
    expected_seeds = ['GGCGGCGGCG', 'CCCGCGCCCC']  # Replace with actual expected sequences if known
    for seed in expected_seeds:
        assert seed in off_target_counter, f"Expected seed sequence '{seed}' not found in off-target hits for Cas3"



def test_extract_sgRNAs_output(sgrna_args_cas9, expected_sgrna_df):
    # Run the extract_sgRNAs function
    sgrna_df = extract_sgRNAs(sgrna_args_cas9)

    # Check that the columns are correct
    expected_columns = expected_sgrna_df.columns
    assert list(sgrna_df.columns) == list(expected_columns), "The columns do not match the expected output"

    # Check the number of rows
    assert len(sgrna_df) == len(expected_sgrna_df), f"The number of rows ({len(sgrna_df)}) does not match the expected ({len(expected_sgrna_df)})"

    # Check that the DataFrame content matches the expected content
    print(sgrna_df)
    print(expected_sgrna_df)
    pd.testing.assert_frame_equal(
        sgrna_df.reset_index(drop=True),
        expected_sgrna_df.reset_index(drop=True),
        atol=1e-5
    )

    # Additional specific checks (optional)
    assert (sgrna_df['gc'] >= sgrna_args_cas9.gc_lower).all() and (sgrna_df['gc'] <= sgrna_args_cas9.gc_upper).all(), "GC content is out of bounds"
    assert sgrna_df['off_target_count'].is_monotonic_increasing, "The DataFrame should be sorted by off-target counts in ascending order"



def test_extract_sgRNAs_cas9_basic(sgrna_args_cas9, expected_sgrna_df):
    # Execute the function
    sgrna_df = extract_sgRNAs(sgrna_args_cas9)

    # Basic checks
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty"
    assert 'strain_name' in sgrna_df.columns, "The resulting DataFrame should contain the strain_name column"
    assert 'gc' in sgrna_df.columns, "The resulting DataFrame should contain the gc column"
    assert 'sgrna' in sgrna_df.columns, "The resulting DataFrame should contain the sgrna column"

    # Check the GC content range
    assert sgrna_df['gc'].between(sgrna_args_cas9.gc_lower, sgrna_args_cas9.gc_upper).all(), "All GC content should be within the specified range"

    # Check off-target sorting
    assert sgrna_df['off_target_count'].is_monotonic_increasing, "The DataFrame should be sorted by off-target counts in ascending order"

    # Compare with the expected DataFrame loaded from CSV
    pd.testing.assert_frame_equal(
        sgrna_df.reset_index(drop=True),
        expected_sgrna_df.reset_index(drop=True),
        atol=1e-5
    )
    print(sgrna_df)
    print(expected_sgrna_df)


def test_extract_sgRNAs_no_hits(sgrna_args_cas9):
    # Modify the arguments to focus on a non-existent gene
    sgrna_args_cas9.locus_tag = ['non_existent_gene']

    # Execute the function
    sgrna_df = extract_sgRNAs(sgrna_args_cas9)

    # The result should be empty since there are no matching genes
    assert sgrna_df.empty, "The resulting DataFrame should be empty for non-existent genes"


def test_extract_sgRNAs_filtering(sgrna_args_cas9):
    # Tighten the GC content filter
    sgrna_args_cas9.gc_upper = 0.5  
    sgrna_args_cas9.gc_lower = 0.4

    # Execute the function with filtering
    sgrna_df = extract_sgRNAs(sgrna_args_cas9)

    # Ensure that all sgRNAs pass the GC content filter
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty"
    assert sgrna_df['gc'].between(0.4, 0.5).all(), "All sgRNAs should have a GC content between 0.4 and 0.5"


def test_extract_sgRNAs_predict_step(sgrna_args_cas9):
    # Add the 'predict' step
    sgrna_args_cas9.step.append('predict')

    # Execute the function
    sgrna_df = extract_sgRNAs(sgrna_args_cas9)

    # Test to make sure the "predict" step runs, though it currently does nothing
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty even if the 'predict' step is included"
    assert 'strain_name' in sgrna_df.columns, "The resulting DataFrame should contain the strain_name column"



def test_extract_sgRNAs_cas3_basic(sgrna_args_cas3, expected_cas3_sgrna_df):
    # Execute the function
    sgrna_df = extract_sgRNAs(sgrna_args_cas3)

    # Basic checks
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty"
    assert 'strain_name' in sgrna_df.columns, "The resulting DataFrame should contain the strain_name column"
    assert 'gc' in sgrna_df.columns, "The resulting DataFrame should contain the gc column"
    assert 'sgrna' in sgrna_df.columns, "The resulting DataFrame should contain the sgrna column"

    # Check the GC content range
    assert sgrna_df['gc'].between(sgrna_args_cas3.gc_lower, sgrna_args_cas3.gc_upper).all(), "All GC content should be within the specified range"

    # Check off-target sorting
    assert sgrna_df['off_target_count'].is_monotonic_increasing, "The DataFrame should be sorted by off-target counts in ascending order"

    # Compare with the expected DataFrame loaded from CSV
    pd.testing.assert_frame_equal(
        sgrna_df.reset_index(drop=True),
        expected_sgrna_df.reset_index(drop=True),
        atol=1e-5
    )
    
    print(sgrna_df)
    print(expected_cas3_sgrna_df)


def test_extract_sgRNAs_no_hits_cas3(sgrna_args_cas3):
    # Modify the arguments to focus on a non-existent gene
    sgrna_args_cas3.locus_tag = ['non_existent_gene']

    # Execute the function
    sgrna_df = extract_sgRNAs(sgrna_args_cas3)

    # The result should be empty since there are no matching genes
    assert sgrna_df.empty, "The resulting DataFrame should be empty for non-existent genes"


def test_extract_sgRNAs_filtering_cas3(sgrna_args_cas3):
    # Tighten the GC content filter
    sgrna_args_cas3.gc_upper = 0.6 
    sgrna_args_cas3.gc_lower = 0.4

    # Execute the function with filtering
    sgrna_df = extract_sgRNAs(sgrna_args_cas3)

    # Ensure that all sgRNAs pass the GC content filter
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty"
    assert sgrna_df['gc'].between(0.4, 0.6).all(), "All sgRNAs should have a GC content between 0.4 and 0.5"


def test_extract_sgRNAs_predict_step_cas3(sgrna_args_cas3):
    # Add the 'predict' step
    sgrna_args_cas3.step.append('predict')

    # Execute the function
    sgrna_df = extract_sgRNAs(sgrna_args_cas3)

    # Test to make sure the "predict" step runs, though it currently does nothing
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty even if the 'predict' step is included"
    assert 'strain_name' in sgrna_df.columns, "The resulting DataFrame should contain the strain_name column"


def test_extract_sgRNAs_output_cas3(sgrna_args_cas3, expected_cas3_sgrna_df):
    # Run the extract_sgRNAs function
    sgrna_df = extract_sgRNAs(sgrna_args_cas3)

    # Check that the columns are correct
    expected_columns = expected_cas3_sgrna_df.columns
    assert list(sgrna_df.columns) == list(expected_columns), "The columns do not match the expected output"

    # Check the number of rows
    assert len(sgrna_df) == len(expected_cas3_sgrna_df), f"The number of rows ({len(sgrna_df)}) does not match the expected ({len(expected_cas3_sgrna_df)})"

    # Check that the DataFrame content matches the expected content
    pd.testing.assert_frame_equal(
        sgrna_df.reset_index(drop=True),
        expected_sgrna_df.reset_index(drop=True),
        atol=1e-5
    )
    print(sgrna_df)
    print(expected_cas3_sgrna_df)

    # Additional specific checks (optional)
    assert (sgrna_df['gc'] >= sgrna_args_cas3.gc_lower).all() and (sgrna_df['gc'] <= sgrna_args_cas3.gc_upper).all(), "GC content is out of bounds"
    assert sgrna_df['off_target_count'].is_monotonic_increasing, "The DataFrame should be sorted by off-target counts in ascending order"


# TODO make a cas12a test suite