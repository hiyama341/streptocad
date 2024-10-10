# Construction and Validation of the Single-Target CRISPR Plasmid

**Timing: 7–16 days**

## 1. Prepare sgRNA Cassettes

1. **Design ssDNA Oligonucleotides**:

   - Design ssDNA oligonucleotides for the single-target CRISPR plasmid construction from the protospacer sequences obtained using the app StreptoCAD. Check the file called : "primers.csv"

2. **Order Oligonucleotides**:

   - Directly order the ssDNA oligonucleotide from the supplier (e.g., IDT) with a standard desalting protocol. You can do this by uploading the csv file called "primers.csv" to IDT.

3. **Resuspend Oligonucleotides**:

   - Resuspend the oligonucleotides in 1× NEBuffer at a 100 μM stock concentration.

4. **Dilute Oligonucleotides**:
   - Before starting the next step, dilute the oligonucleotides to a working concentration of 0.4 μM using 1× NEBuffer.

**PAUSE POINT**: The prepared solutions, together with the stocks, can be stored at −20 °C for up to 6 months.

## 2. Linearize the CRISPR Plasmids by Restriction Enzyme Digestion

**CRITICAL STEP**: Because Thermo Fisher Scientific provides a compatible buffer system, both FastDigest restriction enzymes and FastAP thermosensitive alkaline phosphatase can be added at the same time for a total reaction time of 30 min. Instead of a 20-µL digestion, one can enlarge the volume to yield more linearized plasmid that can be stored at −20 °C for up to 3 months for future use.

### Digesting Non-Multiplexing-Compatible CRISPR Plasmids - ssDNA cloning

1. **Prepare a Digestion**:
   - Prepare a digestion of the selected CRISPR plasmid (pCRISPR–Cas9; pCRISPR–Cas9–ScaligD; pCRISPR–dCas9; pCRISPR–cBEST; or pCRISPR–aBEST) for single-editing applications with FastDigest NcoI in a 20-µL digestion reaction containing the following components:

| Component             | Amount (μL) | Final Concentration |
| --------------------- | ----------- | ------------------- |
| Plasmid DNA           | 10          | 40 ng/µL            |
| FastDigest NcoI       | 1           |                     |
| 10× FastDigest buffer | 2           | 1×                  |
| ddH2O                 | 7           |                     |
| **Total volume**      | 20          |                     |

2. **Digest Plasmid DNA**:
   - Ideally, digest 800 ng of plasmid DNA. Incubate at 37 °C for 30 min. Then add 1 µL of FastAP thermosensitive alkaline phosphatase to the reaction and incubate for an additional 10 min at 37 °C.

**PAUSE POINT**: The linearized plasmids can be stored at −20 °C for up to 3 months.

## 3. Measure the Concentration

Measure the concentration using a NanoDrop 2000 spectrophotometer.

## 4. Insert the sgRNA Cassettes into the Digested CRISPR Plasmid

### Construction of Non-Multiplexing-Compatible CRISPR Plasmids

1. **Prepare Reaction Mix**:
   - To clone the ssDNA oligonucleotide containing the 20-nt protospacer into the linearized single-target CRISPR plasmid using the PCR-free ssDNA oligonucleotide bridging method, first prepare a 20-μL reaction mix as follows:

| Component                                 | Amount (μL) | Final Concentration |
| ----------------------------------------- | ----------- | ------------------- |
| Linearized plasmid                        | 1 (30 ng)   | 1.5 ng/µL           |
| 5 μM oligonucleotide                      | 5           | 1.25 μM             |
| 2x NEBuilder HiFi DNA Assembly Master Mix | 10          | 1×                  |
| ddH2O                                     | 4           |                     |
| **Total volume**                          | 20          |                     |

2. **Incubate Reaction**:
   - Incubate the reaction for 1 h at 50 °C, using a thermocycler.

**CRITICAL STEP**: The reaction volume can be reduced proportionally to 10 μL in order to save reagents.

**PAUSE POINT**: The reaction can be stored at −20 °C for up to 3 months.

3. **Transform E. coli Cells**:
   - Transfer 2 μL of the above reaction mixture to 50 μL of in-house-made (or commercial) electroporation-competent Mach1 E. coli cells. Follow these steps:
     - Remove the 50-µL tubes containing electroporation-competent E. coli cells from the –80 °C freezer and thaw on ice (~10 min).
     - Remove the 1-mm electroporation cuvettes from the −20 °C freezer and place them on ice.
     - Carefully pipette the competent cells into the cuvettes.
     - Mix the competent cells with 2 μL of the ssDNA oligonucleotide bridging reaction by flicking the tubes with a fingertip 3–5 times (avoid bubble formation).
     - Use the Ec1 program of a Bio-Rad MicroPulser (alternatively, a similar electroporation program of 1.8 kV with a 1-mm electroporation cuvette with a one-time pulse can be used).
     - Immediately add 200 µL SOC broth to each cuvette and transfer the reaction to a sterilized 1.5-mL Eppendorf tube.
     - Incubate the tube in a heating block at 37 °C with shaking at 800 r.p.m. for 1 h.

**CRITICAL STEP**: Chemically competent cells could also be used in this step with the 42 °C heat-shock protocol.

4. **Plate Cells**:
   - Plate 100 µL of the reaction onto a selective LB plate supplemented with 50 µg/mL apramycin. Incubate the plate overnight at 37 °C.

**CAUTION**: The transformation efficiency may differ with home-made competent cells. Therefore, we recommend using commercial products if they are available.

## 5. Screen the Clones Using an E. coli Colony PCR

In the morning following the transformation, use sterilized wooden toothpicks to pick 12–24 E. coli colonies from each assembly into a 96-deep-well plate containing 300 µL LB medium supplemented with 50 µg/mL apramycin in each well.

## 6. Incubate the Plate

Incubate the 96-well deep-well plate at 37 °C with 300 r.p.m. shaking for 2 hours.

## 7. Colony PCR Setup and Conditions

Directly use 1 µL of the culture as template DNA for colony PCR with the following setup and conditions:

| Component                                 | Amount (μL)       | Final Concentration |
| ----------------------------------------- | ----------------- | ------------------- |
| sgRNA-TEST-F                              | 0.5               | 500 nM              |
| sgRNA-TEST-R                              | 0.5               | 500 nM              |
| Template DNA                              | 1                 |                     |
| OneTaq 2× Master Mix with standard buffer | 10                | 1×                  |
| ddH2O                                     | 8                 |                     |
| **Total**                                 | 20 (one reaction) |                     |

| Cycle no. | Denature     | Anneal        | Extend       | Final       |
| --------- | ------------ | ------------- | ------------ | ----------- |
| 1         | 94 °C, 3 min |               |              |             |
| 2–31      | 94 °C, 30 s  | \*\* °C, 30 s | 68 °C, 30 s  |             |
| 32        |              |               | 68 °C, 5 min |             |
| 33        |              |               |              | 10 °C, hold |

\*\* Check the annealing temperature from the primers.csv file for the corresponding checking primers.

## 8. Analyze PCR Reaction

Analyze 5 µL of the above PCR reaction (add 1 µL of 6× DNA gel loading dye) along with the GeneRuler 1-kb DNA ladder on a long (10-cm) agarose gel (3% (wt/vol)) with 1× TAE buffer. Run the gel at 100 V for 60 min and visualize the bands using a Gel Doc XR+ Gel Documentation System.

**CRITICAL STEP**: Because the size difference between the positive and the control is only 20 bp, it takes a >2% (wt/vol) agarose gel with a 60-min run time to distinguish the bands. We recommend using a 10-cm 3% (wt/vol) agarose gel and running the gel at 100 V for 60 min.

**PAUSE POINT**: The E. coli culture in the deep-well 96-well plate can be stored at 4 °C for up to 1 week.

## 9. Prepare Cultures of Positive Colonies

Prepare cultures of the above-obtained positive colonies in cultivation tubes containing 5 mL of LB medium supplemented with 50 µg/mL apramycin. Inoculate 50 µL of culture directly from the deep-well 96-well plate and incubate at 37 °C with 200 r.p.m. shaking for ~16 h (overnight).

**CRITICAL STEP**: We always use LB medium for cultivating ET12567 E. coli strains for conjugation but have observed that 2× YT broth performs better than LB medium for plasmid isolation purposes in general. Therefore, we recommend using 5 mL of 2× YT broth for cultivation of each strain that is going to be used for plasmid isolation.

## 10. Perform Plasmid Isolation

Perform plasmid isolation the following day, using the NucleoSpin Plasmid EasyPure Kit and following the manufacturer’s instructions. Submit the plasmids for Sanger sequencing using the sgRNA-TEST-F primer and Cas9-C-terminal-TEST primer. Follow Box 2 to additionally insert the editing template (Fig. 5) after sequencing, if required.

**CRITICAL STEP**: We have observed that instability of the pSG5 replicon–based shuttle plasmid in E. coli is easily triggered by unknown factors. Therefore, it is critical to confirm the integrity of the CRISPR plasmids. We have observed several cases in which the ‘hot region of instability’ lies downstream of the tipA-fd fragment. Therefore, we recommend running an additional Sanger sequencing with the sequencing primer Cas9-C-terminal. Alternatively, a NdeI–BglII double-digestion mapping can also indicate the integrity of the plasmids.

## 11. Freeze the E. coli Strains

Freeze the E. coli strains with the correct plasmids (ready-to-use CRISPR plasmids) confirmed by Sanger sequencing in 25% (vol/vol) glycerol at −80 °C for long-term storage.

**PAUSE POINT**: The E. coli glycerol stock can be stored at −80 °C for at least 5 years.
