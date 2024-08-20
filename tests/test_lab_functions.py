
from streptocad.wet_lab.lab_functions import (
    calculate_master_mix,     
    calculate_volume_and_total_concentration_df, 
    dilute_solution,
    calculate_BsaI_volume
)
import pytest
import pandas as pd
from io import StringIO
import sys
import pandas as pd
from teemi.design.fetch_sequences import read_genbank_files
from pydna.primer import Primer
from pydna.amplify import pcr



def test_calculate_master_mix():
    # Setup
    expected_data = {
        "vol_p_reac": [1.0, 2.5, 2.5, 19.0, 25.0, 50.0],
        "mastermix for 6 reactions": [6.0, "Add primers individually", "Add primers individually", 114.0, 150.0, 300.0]
    }
    expected_df = pd.DataFrame(expected_data, index=["Template", "Primer 1", "Primer 2", "H20", "Pol", "Total"])

    # Call the function
    result_df = calculate_master_mix(
        vol_p_reac=50, 
        no_of_reactions=6, 
        standard_reagents=["Template", "Primer 1", "Primer 2", "H20", "Pol"], 
        standard_volumes=[1, 2.5, 2.5, 19, 25]
    )

    # Assertions
    pd.testing.assert_frame_equal(result_df, expected_df)



def test_calculate_volume_and_total_concentration_df(capsys):
    # Primers
    primers = pd.read_csv('tests/test_files/pcr_df.csv')
    forward_primers = list(primers['f_primer_sequences(5-3)'])
    reverse_primers = list(primers['r_primer_sequences(5-3)'])

    forward_primers_list = [Primer(forward_primers[i]) for i in range(len(forward_primers))]
    reverse_primers_list = [Primer(reverse_primers[i]) for i in range(len(reverse_primers))]
    
    # Template: 
    pJet = read_genbank_files('tests/test_files/Multiplexing_PCR_scaffold.gb')[0]

    # PCR
    amplicon_names = ['PCR_CY00000014', 'PCR_CY00000058', 'PCR_CY00000155', 'PCR_CY00000196', 'PCR_CY00000239']
    concentrations = [55, 72, 25, 78, 103]
    volumes = [45, 45, 45, 45, 45, 45]
    locations = ['l5_D01',  'l5_D02', 'l5_D03',  'l5_D04', 'l5_D05']

    list_of_amplicons = []
    for i in range(len(forward_primers_list)): 
        amplicon = pcr(forward_primers_list[i], reverse_primers_list[i], pJet)
        amplicon.name = amplicon_names[i]
        # Add annotation
        amplicon.annotations['batches'] = [{'location':locations[i], 'concentration':concentrations[i], 'volume': volumes[i]}]
        
        # Save
        list_of_amplicons.append(amplicon)    


    # Prepare input data
    equimolar_concentration_aimed_at = 0.004
    amplicon_names = ('PCR_CY00000014', 'PCR_CY00000058', 'PCR_CY00000155', 'PCR_CY00000196', 'PCR_CY00000239')
    amplicon_parts_amounts_total = {k: v for (k, v) in zip(amplicon_names, [equimolar_concentration_aimed_at] * len(list_of_amplicons))}

    # Expected DataFrame
    expected_df = pd.DataFrame({
        "name": ['PCR_CY00000014', 'PCR_CY00000058', 'PCR_CY00000155', 'PCR_CY00000196', 'PCR_CY00000239'],
        "volume to add": [9.0, 5.7, 16.4, 5.3, 4.0],  # Example values, replace with actual expected values
        "concentration": [55, 72, 25, 78, 103],
        "location": ['l5_D01', 'l5_D02', 'l5_D03', 'l5_D04', 'l5_D05']
    })

    # Call the function
    result_df = calculate_volume_and_total_concentration_df(list_of_amplicons, amplicon_parts_amounts_total)

    print(result_df)
    # Assertions
    pd.testing.assert_frame_equal(result_df, expected_df)
    
    # Capture and assert print output
    captured = capsys.readouterr()

    # Validate printed output
    assert "Total volume of the parts mixed : 40.4" in captured.out
    assert "Final concentration of the parts mixed : 52.990099009900995" in captured.out



def test_dilute_solution():
    # Call the function with provided inputs
    plasmid_to_add, water_to_add = dilute_solution(target_volume=16, concentration=124, final_concentration=10)

    # Assertions to check if the output is as expected
    assert plasmid_to_add == 1.29  # Expecting 1.29 µl of plasmid
    assert water_to_add == 14.71   # Expecting 14.71 µl of water


def test_calculate_BsaI_volume():
    # Test with a typical input where the function should return a valid volume
    # Assuming each fragment adds 0.2 µL to the volume
    number_of_fragments = 5  # This should result in a volume of 1.0 µL
    expected_volume = 1.0
    volume = calculate_BsaI_volume(number_of_fragments)
    assert volume == expected_volume

    # Test with a number of fragments that should raise an error
    # Any number of fragments more than 5 should raise an error, as it would result in more than 1 µL
    with pytest.raises(ValueError) as exc_info:
        calculate_BsaI_volume(6)  # This should result in a volume of 1.2 µL and raise an error

    # Optional: Check if the error message is as expected
    assert "FastDigest BsaI should not exceed 1 µL" in str(exc_info.value)
