# aSPIre
abundance of spliced peptides identified using rescoring

## overview and requirements
*aSPIre* provides the quantification of spliced and non-spliced peptides that were identified with *inSPIRE*.
In order to run *aSPIre*, you need access to a Windows machine that has [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view) installed. **Skyline** is a freely available tool for MS analysis.

*aSPIre* processes peptide-spectrum matches (PSMs) that were assigned by *inSPIRE*, quantifies them using Skyline and constructs a generation kinetic for each identified peptide.


## sample list
The user must provide *inSPIRE* final assignments and features. All information must be provided in the sample_list.csv (see an example below). You can edit the sample list using, e.g., MS Excel. In any case, make sure to save it as file with comma-separated values (.csv) and NOT as .xlsx notebook!

| protein_name | substrateID | substrateSeq | digestTime | biological_replicate | final_assignments | all_features | raw_file |
| ----- |  ----- |  ----- |  ----- |  ----- |  ----- |  ----- |  ----- |
IL37b	| IL37b	| IL37b.fasta| 0 | 	G1	| IL37b_finalAssignments.csv	| IL37b_input_all_features.tab	| WSoh_101121_151121_HFGoe_G1_IL37b_0h_R1.raw
IL37b	| IL37b	| IL37b.fasta	| 4	| G1	| IL37b_finalAssignments.csv	| IL37b_input_all_features.tab	| WSoh_101121_151121_HFGoe_G1_IL37b_4h_R1.raw
IL37b	| IL37b	| IL37b.fasta	| 0	| G2	| IL37b_finalAssignments.csv	| IL37b_input_all_features.tab	| WSoh_101121_151121_HFGoe_G2_IL37b_0h_R1.raw
IL37b	| IL37b	| IL37b.fasta	| 2	| G2	| IL37b_finalAssignments.csv	| IL37b_input_all_features.tab	| WSoh_101121_151121_HFGoe_G2_IL37b_2h_R2.raw

A few general remarks:
- Please always provide a full kinetic including the zero hours / no proteasome control measurements.
- You can put any number of proteins / experiments in the sample list. However, **only one a single can be processed at once**. Specify the protein_name you would like to analyse in `data/config.yaml` (see below).
- Do not put any white spaces, umlauts or empty rows in the sample list. Also, every column has to be filled, even if there are duplicated entries (*e.g.*, substrateID).

### protein_name and substrateID
`protein_name` is the name of the protein/polypeptide to which the respective list entry corresponds. Choose a short and comprehensive name and avoid spaces or special characters.
`substrateID` is the ID under which the given protein will appear in the final output. Can be identical to `protein_name` (recommended), but this does not have to be the case.

### substrateSeq
If you are processing proteins, save their sequence in a single-entry `.fasta` file and deploy it in `data/sequences/`. Specify the name of the `.fasta` file in the `substrateSeq` column of the sample list.

### digestTime
Time point after which the digestion was stopped. **Please provide the time in hours!** For instance, if the digestion time is 15 min, enter 0.25. **Do not put any units in this column!**
Zero-hours time points are treated as control measurements. Please enter *0* as time point for the control measurements.

### biological_replicate
Name of the biological (NOT technical!) replicate. In the final kinetics, the mean over all technical replicates is calculated, whereas biological replicates are displayed separately. Please have a look at the full sample list (`data/sample_list.csv`) for clarification.
Alternatively, instead of the replicate name, you can also put a number in this column. For readability hower, the actual replicate ID is recommended.

### final_assignments and all_features
You will need two files from *inSPIRE* in order to run *aSPIre*: `finalAssignments.csv` and `all_features.tab`. They have to be copied into `data/inSPIRE` and re-named according to the example in the sample list. Specify the name of the *inSPIRE* output files in this column. 

### raw_file
Provide the **full name** (including `.raw` suffix) of the `.raw` file for the respective sample. You do NOT have to copy the `.raw` files to your local device. Make sure that all `.raw` files are accessible to a Windows machine with a Skyline installation.



