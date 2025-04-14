# Conjugations Protocol

**Transfer of Ready-to-Use CRISPR Plasmids into the Target Streptomycetes by Interspecies Conjugation**
**Timing: 2 days**

## 1. Plate the Target Streptomycete Strains

Plate the target streptomycete strains onto MS plates for sporulation. This is typically done days ahead, depending on the growth rate; for example, it requires ~5 days for the S. coelicolor WT strain to fully sporulate on MS plates at 30 °C.

## 2. Transfer Plasmids to E. coli

Transfer 200 ng of the desired plasmids into 50 µL of in-house-made electroporation-competent E. coli ET cells, using the electroporation protocol described previously. Plate onto selective LB plates supplemented with 50 µg/mL apramycin, 50 µg/mL kanamycin, and 25 µg/mL chloramphenicol. Incubate the plates overnight at 37 °C. The next day, wash all transformants off using 4 ml of fresh LB medium, and inoculate liquid LB cultures with 50 ul of cell suspension.

**CRITICAL STEP**: It is also possible to start the culture from a single colony the next day. Using a pool of transformants lowers the chance of picking a wrong clone and can save at least 1 day.

## 3. Prepare ET Cultures

The following morning, prepare the ET cultures harboring the plasmids of interest for conjugation. Harvest 2 ml of overnight culture for each conjugation by centrifugation at 2000g for 2 min, and wash twice using 2 mL of antibiotic-free LB medium.Resuspend the cell pellets in 0.5 mL antibiotic-free LB medium.

## 4. Collect Streptomycete Spores

Collect spores of the streptomycete of interest by pipetting 4 mL 2× YT onto the surface of the spore lawn of a well-sporulated MS plate. Use a sterilized cotton swab to gently scrape off the spores.

## 5. OPTIONAL: Filter Spore Suspension

Carefully place a sterilized cotton pad above the spore suspension with sterilized tweezers.

Use a 25 ml serological pipette to aspirate the spore suspension through the cotton pad to remove agar and mycelial debris; then transfer the suspension to a 50-mL Falcon tube.

**CRITICAL STEP**: Repeat the above operation twice to maximize the amount of spores from one plate. For strains that produce only low numbers of spores, it may be necessary to collect spores from more than one plate.

## 6. OPTIONAL: Heat-Shock Spores

Heat-shock the spore suspension for 10 min at 50 °C. The spore suspension is now ready for conjugation.

**CRITICAL STEP**: To achieve higher conjugation efficiency, leave the spore suspension at 4 °C overnight (up to 3 days) for pre-germination.

**PAUSE POINT**: The spore suspension can be used for up to 2 weeks if stored at 4 °C.

## 7. Mix ET Suspension with Spore Suspension

Mix 500 µL ET suspension with 200 µL spore suspension in a sterilized 1.5-mL Eppendorf tube by carefully pipetting.

## 8. Plate Mixture

Plate the spore-ET cell suspension on MS plates supplemented with 10 mM MgCl2. Air-dry them in a laminar flow hood for 5 min, and incubate the conjugation plates at 30 °C for ~18 h (overnight).

## 9. Overlay Conjugation Plates

Overlay the conjugation plates with 1 mL sterilized H2O containing 250 ug apramycin and 12.5 ug nalidixic acid. Air-dry the plates in a laminar flow hood for 15 min.

**CAUTION**: We do not recommend using a spreader for the overlay procedure. Instead, spread the 1 mL by moving the plate. The 1 mL sterilized H2O will form clouds after you add the nalidixic acid stock.Some streptomycetes can be sensitive to nalidixic acid, in which case it can be left out.

## 10. Incubate Plates

Incubate the plates until exconjugants can be picked with a sterilized wooden toothpick. Typically, it takes 5-7 days if the streptomycete of interest has a normal growth rate. Transfer the picked exconjugants to a fresh ISP2 plate supplemented with 50 µg/mL apramycin and 12.5 µg/mL nalidixic acid and incubate for 5-7 days at 30 °C.

## 11. Evaluate the Successfully Edited Strains

### Evaluation of Non-Genetically Edited Applications (CRISPRi)

1. **Prepare Seed Cultures**:

   - Inoculate three randomly picked exconjugants for each plasmid and a non-treated control into individual shake flasks, each containing 50 mL selective ISP2 broth, and incubate at 30 °C with 180 r.p.m. shaking for 3 days.

2. **Normalize Start Amount**:

   - Spin down 1-mL cultures at 10,000g at room temperature for 5 min to enable measurement of equal amounts of inocula for the main cultures.

3. **Inoculate Main Cultures**:

   - Inoculate 500-mg cell pellets (wet weight) from seed cultures into a fresh 50-mL selective ISP2 broth–containing shake flask and incubate the flask at 30 °C with 180 r.p.m. shaking for 3–5 days.

4. **Induce dCas9**:

   - Induce the cultures with an appropriate amount of thiostrepton. Typically, 0.5 µg/mL is enough to achieve sufficient induction.

5. **Extract and Analyze Endpoint Product**:

   - Extract and analyze the endpoint product accordingly, such as actinorhodin by absorbance measurement.

**CAUTION**: Besides directly analyzing the endpoint product, transcription analysis can also be applied, such as qRT-PCR and RNA-seq.

### Evaluation of Genetically Edited Strains by Colony PCR

1. **Quick Screening**:
   - Scratch ~4 mm2 of mycelia from the ISP2 plates and transfer to a PCR tube containing 40 µL 10 % DMSO.

**CRITICAL STEP**: Use the mycelia from a fast-growing stage (before sporulation) for Streptomyces colony PCR.

2. **Cell Lysis**

   - Microwave the eppendorf tubes with closed lids for 5 min at 300W. If many colonies are lysed at the same time, arrange the tubes in a circle to ensure even microwaving of all samples.
   - Spin down the tubes at max. speed for 5 min.
   - Use 1.5 ul of the supernatant as PCR template.

3. **Colony PCR**:
   - Use 1.5 µL of the supernatant as template DNA for colony PCR. Primers flanking the ~500-bp target regions are used. Perform PCR as follows:

| Component                               | Amount (μL) | Final Concentration |
| --------------------------------------- | ----------- | ------------------- |
| Forward primer                          | 1           | 400 nM              |
| Reverse primer                          | 1           | 400 nM              |
| Template DNA                            | 1.5         | 3% (vol/vol) DMSO   |
| 2× Phusion High-Fidelity PCR Master Mix | 25          | 1×                  |
| ddH2O                                   | 21.5        |                     |
| **Total volume**                        | 50          |                     |

| Cycle no. | Denature     | Anneal      | Extend       | Final       |
| --------- | ------------ | ----------- | ------------ | ----------- |
| 1         | 98 °C, 1 min |             |              |             |
| 2–31      | 98 °C, 10 s  | 65 °C, 30 s | 72 °C, 10 s  |             |
| 32        |              |             | 72 °C, 5 min |             |
| 33        |              |             |              | 10 °C, hold |

**CRITICAL STEP**: We obtained equal efficiencies when using NEB Q5 High-Fidelity 2× Master Mix and 2× Phusion High-Fidelity PCR Master Mix with HF Buffer.

5. **Analyze PCR Products**:

   - Analyze 5 µL of the above PCR reaction along with the GeneRuler 1-kb DNA ladder on an agarose gel (1% (wt/vol)) with 1× TAE buffer. Run the gel at 100 V for 30 min and visualize the bands using a Gel Doc XR+ Gel Documentation System.

6. **Validate Edits by Sanger Sequencing**:
   - If clean PCR bands were obtained, validate the edits by direct Sanger sequencing of 8–12 PCR products, using the forward primers flanking the target regions. Use 0.5 µL of PCR product as template.

### Evaluation of Genetically Edited Mutants by Next Generation Sequencing

We recommend verifying PCR positive mutants by next generation sequencing. For base editing experiments, we recommend using Illumina short read sequencing. For more information, please see:

- https://doi.org/10.1038/s41596-020-0339-z

- https://doi.org/10.1021/acssynbio.3c00188

- https://doi.org/10.1073/pnas.1913493116

For any structural mutations such as deletions and integrations, we recommend long read sequencing such as Oxford Nanopore sequencing. For more information, please see:

- https://doi.org/10.1016/j.xpro.2022.101955

- https://doi.org/10.1093/nar/gkaf214
