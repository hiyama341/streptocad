# Multi-Target CRISPR Protocol

**Timing: 7–16 days**

### Preparing sgRNA Cassettes for Multiple-Target CRISPR Plasmid Construction

1. **sgRNA Cassette Preparation**:
   - sgRNA cassette preparation for multiplexing-compatible CRISPR plasmids requires PCR to obtain the cassettes. First, design and order primers for each sgRNA fragment according StreptoCAD file generate called "primers.csv".

**CAUTION**: Check that the overhang dataframe does not contain duplicates.

2. **PCR for Each sgRNA**:
   - Carry out a PCR for each sgRNA in a 50-µL reaction, using the pJET1.2–sgRNA handle as template, with the following setup and conditions:

| Component                               | Amount (μL)  | Final Concentration |
| --------------------------------------- | ------------ | ------------------- |
| Forward primer                          | 1            | 400 nM              |
| Reverse primer                          | 1            | 400 nM              |
| Template DNA                            | 0.5 (~50 ng) |                     |
| 2× Phusion High-Fidelity PCR Master Mix | 25           | 1×                  |
| ddH2O                                   | 22.5         |                     |
| **Total**                               | 50           |                     |

| Cycle no. | Denature    | Anneal        | Extend       | Final       |
| --------- | ----------- | ------------- | ------------ | ----------- |
| 1         | 98 °C, 20 s |               |              |             |
| 2–31      | 98 °C, 10 s | 60\* °C, 30 s | 72 °C, 10 s  |             |
| 32        |             |               | 72 °C, 5 min |             |
| 33        |             |               |              | 10 °C, hold |

**CAUTION**: The extension time is calculated according to the DNA polymerase used, for example, 15–30 s/kb for Phusion High-Fidelity polymerase. \* Use the correct value from the "primers.csv" file.

**CRITICAL STEP**: We obtained equal efficiencies by using the NEB Q5 High-Fidelity 2× Master Mix Kit and the 2× Phusion High-Fidelity PCR Master Mix with HF Buffer.

**CRITICAL STEP**: The annealing temperature (Ta) is calculated with a Tm calculator from NEB (https://tmcalculator.neb.com/#!/main).

3. **Analyze PCR Reaction**:
   - Analyze 5 µL of each PCR reaction (add 1 µL of 6× DNA gel loading dye) along with the GeneRuler 100-bp DNA ladder on an agarose gel (1% (wt/vol)) with 1× TAE buffer. Run the gel at 100 V for 30 min and visualize the bands using a Gel Doc XR+ Gel Documentation System. A band at ~150 bp is expected.

**TROUBLESHOOTING**

4. **Purify Positive Fragments**:
   - Purify the positive fragments using a NucleoSpin Gel and PCR Clean-up Kit according to the manufacturer’s instructions.

**CRITICAL STEP**: The reaction can be directly used for the following digestion–ligation processes; however, we recommend first purifying the fragments.

5. **Measure Fragment Concentration**:

   - Measure the concentration of each fragment using a NanoDrop 2000 spectrophotometer.

6. **Pool Fragments**:
   - Pool all fragments, calculate the volume of each fragment required based on the measured concentrations and a desired molar ratio of 5:1 (insert/vector).

**CRITICAL STEP**: This equation is used for picomole-to-nanogram conversion: Picomoles = (weight in nanograms) × 1,000/(base pairs × 650 Da).

7. **Digest sgRNA Fragments**:
   - Digest the sgRNA fragments using FastDigest Eco31I (BsaI). Set up a 20-µL reaction using a mix of 1 µL FastDigest Eco31I and 1 µL FastAP thermosensitive alkaline phosphatase and incubate for 30 min at 37 °C. Heat-inactivate for 10 min at 65 °C.

**CRITICAL STEP**: The amount of FastDigest Eco31I (BsaI) needs to be calculated based on the number of fragments; we recommend using 0.2 µL FastDigest Eco31I for each fragment in a 20-µL reaction. In total, the amount of FastDigest Eco31I should not exceed 1 µL. If there are more than 5 fragments, use a NEB BsaI-HF v2 kit to insert the sgRNA array into the pGGA vector included in the kit and then apply a NcoI–NheI double digestion to isolate the pre-assembled sgRNA array.

**CRITICAL STEP**: Because there are many BsaI sites in the pSG5 plasmid backbone, a two-step Golden Gate assembly reaction must be set up to insert the sgRNA fragments into the multiplexing-compatible CRISPR plasmids. Therefore, removal of all restriction enzymes by proper heat inactivation or cleanup of the fragments is important.

**PAUSE POINT**: The prepared sgRNA fragments can be stored at −20 °C for up to 6 months.

## 2. Digest Multiplexing-Compatible CRISPR Plasmids

1. **Double Digestion**:
   - Use NcoI and NheI to digest the multiplexing-compatible CRISPR plasmids. Prepare a double digestion of the multiplexing-compatible CRISPR plasmid (pCRISPR–McBEST) with FastDigest NcoI and FastDigest NheI in a 20-µL digestion reaction containing the following components:

| Component             | Amount (μL) | Final Concentration |
| --------------------- | ----------- | ------------------- |
| Plasmid DNA           | 10          | 40 ng/µL            |
| FastDigest NcoI       | 1           |                     |
| FastDigest NheI       | 1           |                     |
| 10× FastDigest buffer | 2           | 1×                  |
| ddH2O                 | 6           |                     |
| **Total volume**      | 20          |                     |

2. **Digest Plasmid DNA**:
   - Ideally, digest 800 ng of plasmid DNA. Incubate it at 37 °C for 30 min and then add 1 µL of FastAP thermosensitive alkaline phosphatase to the reaction and incubate for an additional 10 min at 37 °C, followed by an inactivation step at 75 °C for 10 min.

**PAUSE POINT**: The linearized plasmids can be stored at −20 °C for up to 3 months.

## 3. Analyze Digestion Reaction

Analyze 2 µL of the above digestion reaction (add 1 µL of 6× DNA gel loading dye and 4 µL of ddH2O) along with the GeneRuler 1-kb DNA ladder on an agarose gel (1% (wt/vol)) with 1× TAE buffer. Run the gel at 100 V for 30 min and visualize the bands using a Gel Doc XR+ Gel Documentation System.

## 4. Clean Up Linearized Plasmid

Clean up the gel-confirmed linearized plasmid with a NucleoSpin Gel and PCR Clean-up Kit according to the manufacturer’s instructions.

**CRITICAL STEP**: Optionally, a gel purification step can also be used for extraction of the linearized plasmid; however, it will typically give a much lower yield.

## 5. Measure the Concentration

Measure the concentration using a NanoDrop 2000 spectrophotometer.

## 6. Insert sgRNA Cassettes into the Multiplexing-Compatible CRISPR Plasmids by a Two-Step Golden Gate Assembly

1. **Set Up Ligation Reaction**:
   - Set up a 10-µL ligation reaction using T4 ligase (5 U), 100 ng of predigested plasmids, and the required volume of the prepared sgRNA fragments to give a 5:1 molar ratio of insert/vector (use the heat-inactivated digestion reaction). Supplement the reaction with 2 µL of 50% (wt/vol) PEG-4000 (included in the T4 ligase kit) and incubate for 1 h at 22 °C.

**PAUSE POINT**: The reaction can be stored at −20 °C for up to 3 months.

**TROUBLESHOOTING**

2. **Transform E. coli Cells**:
   - Transfer 5 µL of the above reaction to 50 µL One Shot Mach1 T1 phage-resistant chemically competent E. coli, using a heat-shock protocol as follows: remove 50-µL aliquots of the chemically competent E. coli cells from the –80 °C freezer and thaw on ice (~10 min). Mix with 5 µL of the reaction by flicking with a fingertip (3–5 times). Incubate on ice for 20 min, followed by a 60-s heat shock at 42 °C, using a water bath; transfer the tubes into the ice and incubate for another 5 min. Add 200 µL SOC to the tubes, incubate in a heating block under conditions of 37 °C with 800 r.p.m. shaking for 1 h.

**CRITICAL STEP**: Electroporation-competent cells could also be used in this step with an electroporation transformation protocol.

3. **Plate Cells**:
   - Plate all of the reaction on a selective LB plate supplemented with 50 µg/mL apramycin. Incubate the plate overnight at 37 °C.

## 7. Screen Clones Using E. coli Colony PCR

In the morning following the transformation, use sterilized wooden toothpicks to pick 12–24 E. coli colonies from each assembly into a 96-deep-well plate containing 300 µL LB medium supplemented with 50 µg/mL apramycin in each well.

## 8. Incubate Plate

Incubate the 96-well deep-well plate at 37 °C with 300 r.p.m. shaking for 2 hours.

## 9. Colony PCR Setup and Conditions

Directly use 1 µL of the culture as template DNA for colony PCR with the following setup and conditions:

| Component                                 | Amount (μL)       | Final Concentration |
| ----------------------------------------- | ----------------- | ------------------- |
| sgRNA-TEST-F                              | 0.5               | 500 nM              |
| sgRNA-TEST-R                              | 0.5               | 500 nM              |
| Template DNA                              | 1                 |                     |
| OneTaq 2× Master Mix with standard buffer | 10                | 1×                  |
| ddH2O                                     | 8                 |                     |
| **Total**                                 | 20 (one reaction) |                     |

| Cycle no. | Denature     | Anneal      | Extend       | Final       |
| --------- | ------------ | ----------- | ------------ | ----------- |
| 1         | 94 °C, 3 min |             |              |             |
| 2–31      | 94 °C, 30 s  | 50 °C, 30 s | 68 °C, 30 s  |             |
| 32        |              |             | 68 °C, 5 min |             |
| 33        |              |             |              | 10 °C, hold |

## 10. Analyze PCR Reaction

Analyze 5 µL of the above PCR reaction (add 1 µL of 6× DNA gel loading dye) along with the GeneRuler 1-kb DNA ladder on a long (10-cm) agarose gel (3% (wt/vol)) with 1× TAE buffer. Run the gel at 100 V for 60 min and visualize the bands using a Gel Doc XR+ Gel Documentation System.

**CRITICAL STEP**: Because the size difference between the positive and the control is only 20 bp, it takes a >2% (wt/vol) agarose gel with a 60-min run time to distinguish the bands. We recommend using a 10-cm 3% (wt/vol) agarose gel and running the gel at 100 V for 60 min.

**PAUSE POINT**: The E. coli culture in the deep-well 96-well plate can be stored at 4 °C for up to 1 week.

## 11. Prepare Cultures of Positive Colonies

Prepare cultures of the above-obtained positive colonies in cultivation tubes containing 5 mL of LB medium supplemented with 50 µg/mL apramycin. Inoculate 50 µL of culture directly from the deep-well 96-well plate and incubate at 37 °C with 200 r.p.m. shaking for ~16 h (overnight).

**CRITICAL STEP**: We always use LB medium for cultivating ET12567 E. coli strains for conjugation but have observed that 2× YT broth performs better than LB medium for plasmid isolation purposes in general. Therefore, we recommend using 5 mL of 2× YT broth for cultivation of each strain that is going to be used for plasmid isolation.

## 12. Perform Plasmid Isolation

Perform plasmid isolation the following day, using the NucleoSpin Plasmid EasyPure Kit and following the manufacturer’s instructions. Submit the plasmids for Sanger sequencing using the sgRNA-TEST-F primer and Cas9-C-terminal-TEST primer. Optionally, follow additional steps to insert the editing template after sequencing, if required.

**CRITICAL STEP**: We have observed that instability of the pSG5 replicon–based shuttle plasmid in E. coli is easily triggered by unknown factors. Therefore, it is critical to confirm the integrity of the CRISPR plasmids. We have observed several cases in which the ‘hot region of instability’ lies downstream of the tipA-fd fragment. We recommend running an additional Sanger sequencing with the sequencing primer Cas9-C-terminal. Alternatively, a NdeI–BglII double-digestion mapping can also indicate the integrity of the plasmids.

## 13. Freeze E. coli Strains

Freeze the E. coli strains with the correct plasmids (ready-to-use CRISPR plasmids) confirmed by Sanger sequencing in 25% (vol/vol) glycerol at −80 °C for long-term storage.

**PAUSE POINT**: The E. coli glycerol stock can be stored at −80 °C for at least 5 years.
