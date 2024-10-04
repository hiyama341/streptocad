# Integration and Validation of genes of interest into the overexpression pOEX_PkasO Plasmid

**Timing: 14–21 days**

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

   - Inoculate E. coli ET122567 from a cryotube into 5 mL of YT medium supplemented with 50 µg/mL Kanamycin and 25 µg/mL Chloramphenicol. Incubate overnight at 37°C with shaking.

2. **Day 2: Scale-up Culture**
   - Inoculate a 500 mL flask containing 2×YT medium with the pre-culture to reach an OD600 of approximately 0.075.
   - Supplement with Kanamycin (50 µg/mL) and Chloramphenicol (25 µg/mL).
   - Incubate at 37°C with shaking until the culture reaches an OD600 of 0.45–0.5.

## Harvest and Wash Cells for Electroporation

1. **Harvest Cells:**

   - Aliquot 1.4 mL of the culture into microcentrifuge tubes and centrifuge at 4000×g for 1 min.
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

**PAUSE POINT**: Electroporated cells can be stored at 4°C for up to 24 hours before plating.

## Plate Cells and Select Transformants

1. **Plate Cells:**

   - Plate 100 µL of the electroporated cells onto LB agar supplemented with Kanamycin (50 µg/mL), Chloramphenicol (25 µg/mL), and Apramycin (12.5 µg/mL).
   - Incubate overnight at 37°C.

2. **Pick Colonies:**
   - Pick 6–12 colonies from each transformation and inoculate in 5 mL of LB with the same antibiotics. Incubate overnight at 37°C with shaking.

## Colony PCR to Verify Integration

1. **Prepare Colonies:**

   - Transfer a small amount of each colony into 50 µL of 10% DMSO for lysis.

2. **Cell Lysis:**

   - Freeze cells on dry ice for 10 min, microwave for 2 min, and incubate at 99°C for 10 min. Repeat the freezing and heating cycle once more.

3. **Colony PCR Setup:**

| Component | Volume (µL) |
| --------- | ----------- |
| Q5 (GC)   | 5           |
| Primer F  | 0.5         |
| Primer R  | 0.5         |
| ddH₂O     | 4           |
| **Total** | 10          |

- **PCR Cycling Conditions:**

| Step  | Temperature | Time      |
| ----- | ----------- | --------- |
| 1     | 98°C        | 30 sec    |
| 2–30  | 98°C        | 10 sec    |
|       | 64°C        | 30 sec    |
|       | 72°C        | 60 sec/kb |
| Final | 72°C        | 5 min     |

4. **Analyze PCR Products:**
   - Run 5 µL of PCR products on a 1% agarose gel.

**PAUSE POINT**: PCR products can be stored at −20°C for up to 1 week.

## Next step is conjugation. Please refer to the conjugation_protocol
