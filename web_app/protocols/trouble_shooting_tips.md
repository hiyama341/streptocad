# Streptomyces Protocols Troubleshooting Guide

This guide is based on the tips from the paper: https://www.nature.com/articles/s41596-020-0339-z#Sec60. If you need help with anything refer to this paper or write an email :)

## 1. PCR Amplification

### Problem: No Bands

**Cause:** Forgot to add primers or DNA template  
**Solution:** Double-check that both primers and template are present in the mix; repeat PCR

### Problem: Unspecific or Multiple Bands

**Cause:** High-GC templates and long primer overhangs can promote off-target binding  
**Solutions:**

- Add 3% (vol/vol) DMSO to the reaction
- Use a **touchdown PCR** program (we have found it effective to start 6°C above the calculated annealing temperature, then decrease three times by 2°C, each cycle). For more info: https://www.neb.com/en/faqs/0001/01/01/what-is-touchdown-pcr

---

## 2. Colony PCR

### Problem: No Amplification

**Cause:** Insufficient DNA released from colony lysate  
**Solutions:**

- Pick slightly more cells (but avoid carryover of agar)
- Alternatively, use 1–10 ng of purified genomic DNA as template

### Problem: Strong Background or Smears

**Cause:** Debris in crude lysate inhibits specificity  
**Solutions:**

- Reduce lysate volume (e.g., 1 µL instead of 2–3 µL)
- Switch to a high-GC-optimized polymerase or apply touchdown PCR

---

## 3. DNA Ligation/Hifi-Gibson assembly

### Problem: No Colonies After Transformation

**Cause:** Impure fragments or incompatible overhangs  
**Solutions:**

- **Gel-purify** both insert and vector
- Confirm the ratios are > 1.8
- Verify overhang compatibility using an in silico tool such as Pydna.

---

## 4. Conjugation into Streptomyces

### Problem: Too Few or No Exconjugants

**Causes:**

- Low spore or donor cell density
- Plasmid instability or incompatible replicon
- Media composition issues

**Solutions:**

- Increase both spore count and _E. coli_ donor volume
- Sequence the pSG5 (or alternative) replicon region from the donor to ensure integrity
- If pSG5 is incompatible, swap in a compatible Streptomyces replicon
- Prepare media carefully—avoid full-fat soy flour and add MgCl₂ **after** autoclaving

---

## 5. CRISPRi (dCas9) Knockdown

### Problem: Little or No Knockdown

**Cause:** Ineffective guide RNA

**Solutions:**

- Design and test at least three different protospacers
- Replace the tipA promoter with a known strong promoter if necessary
- Measure transcript reduction by qRT-PCR or RNA-seq at an appropriate time point rather than relying on endpoint phenotypes

---

## 6. Sequencing Validation

### Problem: Mixed or Noisy Sequencing Peaks

**Cause:** Mixed population of edited and unedited cells

**Solutions:**

- Re-streak single colonies on selective plates
- Re-isolate DNA and re-sequence
- If unwanted signal is < 20%, the edit may still be functional—but verify by phenotype and sequencing

### Problem: Low Editing Efficiency

**Cause:** Suboptimal spacer or low Cas9/dCas9 expression

**Solutions:**

- For inducible systems, include inducer (e.g., thiostrepton) during conjugation and outgrowth
- Try alternative protospacers if efficiency remains low

---

## 7. Genomic DNA Isolation

### Problem: Low or No DNA Yield

**Cause:** Overloaded extraction columns  
**Solution:** Split biomass into two preparations or halve the starting material

### Problem: Incomplete Lysis

**Cause:** Insufficient lysozyme or incubation time  
**Solution:** Double lysozyme concentration and extend lysis incubation

### Problem: Degraded DNA

**Cause:** DNase contamination or harsh handling  
**Solutions:**

- Use fresh reagents
- Avoid vortexing; pipette gently to minimize shearing

---

## Quick Reference

| Protocol        | Most Common Issue | First Step to Try                    |
| --------------- | ----------------- | ------------------------------------ |
| PCR             | No bands          | Check primers and template are added |
| Colony PCR      | No amplification  | Pick more cells, avoid agar          |
| Ligation/Gibson | No colonies       | Gel-purify insert and vector         |
| Conjugation     | Few exconjugants  | Increase cell densities              |
| CRISPRi         | Poor knockdown    | Test multiple protospacers           |
| Sequencing      | Mixed peaks       | Re-streak single colonies            |
| gDNA isolation  | Low yield         | Reduce starting material             |
