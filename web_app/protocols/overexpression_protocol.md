# Integration and Validation of genes of interest into the overexpression pOEX_PkasO Plasmid

**Timing: 10-14 days**

## Amplification of Genes of Interest

### 1. Amplify Genes of Interest Using PCR

1. **Prepare the PCR Reaction:**

   - Use the `primer_df` file as a guide to identify the correct primer pairs for each gene of interest. Each row in the `primer_df` corresponds to a specific target gene with the necessary forward and reverse primer sequences.

2. **PCR Reaction Setup:**

| Component                       | Volume (µL) | Final Concentration |
| ------------------------------- | ----------- | ------------------- |
| 5× Q5 Reaction Buffer           | 10          | 1×                  |
| dNTP Mix (10 mM)                | 1           | 200 µM              |
| Forward Primer (10 µM)          | 2.5         | 0.5 µM              |
| Reverse Primer (10 µM)          | 2.5         | 0.5 µM              |
| Template DNA                    | 1           | ~1-50 ng            |
| Q5 High-Fidelity DNA Polymerase | 0.5         | -                   |
| Nuclease-free water             | Up to 25 µL |                     |

3. **PCR Cycling Conditions:**

| Step                 | Temperature          | Time      | Cycles               |
| -------------------- | -------------------- | --------- | -------------------- |
| Initial Denaturation | 98°C                 | 30 sec    | 1                    |
| Denaturation         | 98°C                 | 10 sec    | 30                   |
| Annealing            | Refer to `primer_df` | 30 sec    | Refer to `primer_df` |
| Extension            | 72°C                 | 30 sec/kb | 30                   |
| Final Extension      | 72°C                 | 5 min     | 1                    |

4. **Verify PCR Products:**
   - Run 5 µL of each PCR product on a 1% agarose gel alongside a DNA ladder to confirm the expected size of the amplified genes.

**PAUSE POINT**: Amplified PCR products can be stored at −20°C for up to 1 week.

## Gel extractions##

1. If the band size is correct, run a gel for gel extractions.
2. Purify the correct band from the gel using your preferred kit. For Macherey-Nagel kits, incubate the colum with elution buffer in a thermoblock for 5 min at 70 C prior to elution to increase the recovery.
3. Measure the fragment concentration on a Nanodrop.

## Plasmid Assembly##

1. **Digestion of pOEX-kasOP**
   - Digest 1 ug of pOEX-kasOP using StuI and FastAP.
   - Incubate for 2 h at 37 C, followed by inactivation for 10 min at 75 C.

| Component             | Volume (µL) | Final Concentration |
| --------------------- | ----------- | ------------------- |
| pEOX-kasOP            | Calculate   | 1 ug                |
| 10x FastDigest buffer | 6           | 1x                  |
| StuI                  | 2           | 20 U                |
| FastAP                | 2           | 2 U                 |
| ddH2O                 | to 60 ul    | -                   |

3. **Gibson Assembly**
   - Set up a Gibson Assembly using the digested pOEX-kasOP and the gel extracted fragment.
   - Incubate at 50 C for 1h.

| Component                      | Volume (µL) | Final Concentration |
| ------------------------------ | ----------- | ------------------- |
| pEOX-kasOP*StuI*FastAP         | Calculate   | 100 ng              |
| Gene Fragment                  | Calculate   | 3x equimolar mass   |
| 2x Hifi DNA Assembly Mastermix | 5           | 1x                  |
| ddH2O                          | to 10 ul    | -                   |

4. **Transformations into E. coli Mach1**

   - Transform 4 ul of the assembly mix into E. coli Mach1 or equivalent cloning strains using chemical transformation, following the manufacturers protocol.
   - Plate the transformation on LB plates supplemented with 50 ug/ml Apramycin.
   - Incubate at 37 C overnight.

5. **Colony PCRs**

   - Pick transformants on a fresh LB+Apr plate using wooden toothpicks. After picking, add the toothpicks to PCR tube with 20 ul of ddH20. Twist the toothpicks to dispense remaining biomass into the ddH2O.
   - Use 1 ul of the ddH2O cell suspension as template for colony PCRs.

     PCR Protocol:
     Now we need to check if the gene of interest is correctly integrated into the plasmid. We will do this by PCR. We have a forward and reverse primer for this.

   - f_primer : gcggtgttgtaaagtcgtggcc
   - r_primer : ccgatcaaccgcgactagcatcg

   Since we use the same primers for all the pcrs we can group the amplicons according to size (keep temp really high).

   We typically use a touchdown PCR protocol with three distinct steps:

   1. 72 C for x sec according to the group elongation time
   2. 70 C for x sec according to the group elongation time
   3. 66 C for x sec according to the group elongation time

   Note: If you have a lot of PCR to perform it can be advantagous to group amplicons with similar sizes and in the same thermocycler run. We use the function group_amplicons() to group the amplicons into a dataframe with amplicons with similar sizes (See for inspiration: https://github.com/hiyama341/streptocad/blob/main/notebooks/wet_lab_notebooks/02-Integration_of_G%C3%964010_regulators_into_pOEX_PkasO.ipynb)

6. **Overnight cultures**
   - Prepare overnight cultures of the positive colonies using the plate with picked colonies. Use culture tubes. Inoculate 4 ml of 2x YT supplemented with 50 ug/ml of apramycin directly from the plate.
   - Incubate overnight at 37 C while shaking.
7. **Minipreps and glycerol stocks**
   - Perform minipreps using 2 ml of overnight culture using your preferred kit. Follow the manufacturers instructions.
   - Keep the remaining culture in the cold room until the next day to prepare glycerol stocks of the sequence verified clones.
   - Glycerol stocks are prepared using equal volume of culture and 50 % sterile glycerol.
8. **Sanger sequencing**
   - For sequence verification, submit the miniprepped plasmid using your preferred sequencing kit. Follow the manufacturers instructions. Submit 2 samples for each plasmid, one with CW1026 and one with CW1027.

## Pre-culture Preparation of E. coli ET122567

1. **Day 1: Start Pre-culture**

   - Inoculate E. coli ET12567 from a cryotube into 5 mL of YT medium supplemented with 50 µg/mL Kanamycin and 25 µg/mL Chloramphenicol. Incubate overnight at 37°C with shaking.

2. **Day 2: Scale-up Culture**
   - Inoculate a 500 mL flask containing 2×YT medium with the pre-culture to reach an OD600 of approximately 0.075.
   - Supplement with Kanamycin (50 µg/mL) and Chloramphenicol (25 µg/mL).
   - Incubate at 37°C with shaking until the culture reaches an OD600 of 0.45–0.5.

## Harvest and Wash Cells for Electroporation

1. **Harvest Cells:**

   - Aliquot 2 mL of the culture for each transformation into microcentrifuge tubes and centrifuge at 4000×g for 1 min.
   - Discard the supernatant and wash the cell pellet twice with 500 µL of ddH₂O, centrifuging at the same speed each time.

2. **Resuspend Cells:**
   - Resuspend each cell pellet in 50 µL of sterile ddH₂O.

## Electroporation of E. coli with pOEX_PkasO

1. **Mix DNA with Competent Cells:**

   - Add 1–2 µL of plasmid DNA to the resuspended cells and gently mix.
   - Transfer the mixture into a 1-mm electroporation cuvette.

2. **Perform Electroporation:**

   - Use the EC1 program on a Bio-Rad MicroPulser to perform electroporation.
   - Immediately add 1 mL of 2×YT medium and transfer the cells to a sterile tube.

3. **Incubate:**
   - Incubate the cells for 1–1.5 hours at 37°C with shaking.

## Plate Cells and Harvest Transformants

1. **Plate Cells:**

   -Spin down the cells at 4000xg for 2 min and decant the supernatant. Resuspend the cells in the flowback.

   - Plate the electroporated cells onto LB agar supplemented with Kanamycin (50 µg/mL), Chloramphenicol (25 µg/mL), and Apramycin (12.5 µg/mL).
   - Incubate overnight at 37°C.

2. **Cell Harvesting:**
   - Wash all transformants off the plate using an L-spreader and 2-4 ml of 2xYT medium. Use some of the cell suspension to prepare a glycerol stock, and use the left over to inoculate an overnight culture in 5 ml of 2xYT supplemented with Kanamycin (50 µg/mL), Chloramphenicol (25 µg/mL), and Apramycin (12.5 µg/mL). Incubate overnight at 37°C with shaking.

## Conjugations##

To be performed in a LAF bench.

1. **ET cell preparation**

   - Harvest 2 ml of the ET + plasmid overnight cultures and spin down at 2000xg for 2 min. Gently remove the supernatant and wash once with 1ml of 2xYT medium.
   - Resuspend in 500 ul of 2xYT medium

2. **Mixing Spores and ETs**

   - Mix 500 ul of resuspended ET+plasmid cells with 500 ul of spores of the Streptomyces strain of interest.
   - Plate on the appropriate medium for conjugation (default for us: Mannitol Soy Flour medium supplemented with 10 mM MgCl2)
   - Let the plates dry in the airflow of the LAF bench. Once dried, incubate facing upwards at 30 C.

3. **Overlay with Antibiotics**
   - After 18-24 h, overlay the plates with 1ml of ddH2O supplemented with 5 ul of 50 mg/ml of Apramycin. Spread without using a spreader by just carefully tilting the plate.
   - Let the plate dry fully in the airflow of the LAF bench, and rotate occasionally.
   - Once fully dried, return the plates into the incubator and continue incubation at 30 C. Turn the plates upside down now.
4. **Colony Picking**
   - After 7-10 days (depening on the strain), pick exconjugants using sterile wooden toothpicks and transfer to ISP2 plates supplemented with 50 ug/ml of Apramycin and 25 ug/ml Nalidixic Acid.
   - Incubate the plates until fully grown at 30 C.

## Colony PCR to Verify Integration

1. **Prepare Colonies:**

   - Transfer a small amount of each colony into 50 µL of 10% DMSO for lysis in an eppendorf tube.

2. **Cell Lysis:**

   - Microwave the eppendorf tubes with closed lids for 5 min at 300W. If many colonies are lysed at the same time, arrange the tubes in a circle to ensure even microwaving of all samples.
   - Spin down the tubes at max. speed for 5 min.
   - Use 1.5 ul of the supernatant as PCR template.

3. **Colony PCR Setup:**

| Component | Volume (µL) |
| --------- | ----------- |
| Q5 (GC)   | 5           |
| Primer F  | 0.5         |
| Primer R  | 0.5         |
| ddH₂O     | 4           |
| **Total** | 10          |

- **PCR Cycling Conditions:**

| Step    | Temperature        | Time      |
| ------- | ------------------ | --------- |
| Initial | 98°C               | 30 sec    |
| 2-30    | 98°C               | 10 sec    |
|         | Refer to primer_df | 30 sec    |
|         | 72°C               | 60 sec/kb |
| Final   | 72°C               | 5 min     |

4. **Analyze PCR Products:**
   - Run 5 µL of PCR products on a 1% agarose gel.

**PAUSE POINT**: PCR products can be stored at −20°C for up to 1 week.

## Next step is conjugation. Please refer to the conjugation_protocol
