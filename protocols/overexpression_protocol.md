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

4. **Overlay with Antibiotics**
   - After 18-24 h, overlay the plates with 1ml of ddH2O supplemented with 5 ul of 50 mg/ml of Apramycin. Spread without using a spreader by just carefully tilting the plate.
   - Let the plate dry fully in the airflow of the LAF bench, and rotate occasionally.
   - Once fully dried, return the plates into the incubator and continue incubation at 30 C. Turn the plates upside down now.
     
3. **Colony Picking**
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
