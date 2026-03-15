# Construction and Validation of the Single-Target CRISPR-Cas3 Plasmid

**Timing: 7-16 days**

## 1. Order the Cas3 Spacer-Cloning Primers

1. Use the numbered IDT order file from StreptoCAD (`01_full_idt.csv` or `02_full_idt.csv`, depending on whether repair-template primers are included).
2. Order the spacer-specific primers listed for your selected protospacer(s).
3. Also order the two fixed Cas3 backbone cloning primers:
   - `Cas3_backbone_cloning_fwd`: `GAGCTCATAAGTTCCTATTCCGAAG`
   - `Cas3_backbone_cloning_rev`: `AAGAAGTGGGTGTCGGACGC`

The spacer-specific primers carry the 34-nt protospacer in their overhangs, and the fixed backbone primers amplify the two plasmid fragments used for Gibson assembly.

## 2. Linearize the pCRISPR-Cas3 Backbone

Digest the pCRISPR-Cas3 backbone with `BstBI` and `NdeI`, then purify the digested backbone before Gibson assembly.

## 3. Amplify the Two Spacer-Containing PCR Fragments

For each protospacer, perform two PCRs on the pCRISPR-Cas3 backbone:

- PCR 1: fixed backbone forward primer + spacer-specific reverse primer
- PCR 2: spacer-specific forward primer + fixed backbone reverse primer

StreptoCAD already designs the spacer-specific primers for this Gibson strategy. The two fixed backbone primers are also included in the numbered IDT order file (`01_full_idt.csv` or `02_full_idt.csv`).

## 4. Assemble the Spacer into pCRISPR-Cas3 by Gibson Assembly

Set up a Gibson assembly with:

- the `BstBI + NdeI` digested pCRISPR-Cas3 backbone
- PCR fragment 1
- PCR fragment 2

Incubate using your standard NEBuilder HiFi / Gibson conditions.

## 5. Transform and Screen E. coli Clones

Transform the Gibson assembly into a cloning strain such as Mach1 or an equivalent E. coli strain and plate on apramycin selection.

For colony PCR screening of spacer-positive clones, use these fixed screening primers:

- `CW1-sgRNA-seq_fwd`: `GTACGCGGTCGATCTTGACG`
- `CW2-sgRNA-seq_rev`: `TGCTGACCGGATCAGCAGTC`

Sequence positive clones to confirm the spacer insert and plasmid integrity.

## 6. Clone Repair Templates in a Second Round if Needed

If you are building an in-frame deletion construct, clone repair templates into the MCS after the spacer cloning step.

StreptoCAD provides a workflow-order file to help determine cloning order. This matters if your repair template contains `BstBI` or `NdeI` sites, in which case the repair template cloning should be done after spacer cloning.

## 7. Recommended Interpretation of StreptoCAD Outputs

- `01_full_idt.csv` or `02_full_idt.csv`: order list for spacer primers, fixed Cas3 backbone cloning primers, and checking primers
- `filtered_sgrna_df.csv`: selected protospacers
- `plasmid_metadata_df.csv`: plasmid overview
- `workflow_order_df.csv`: recommended cloning order when repair templates are included

For workflow 6, follow this Cas3-specific Gibson cloning protocol rather than the ssDNA-bridging protocol used for the single-target Cas9-style plasmids.
