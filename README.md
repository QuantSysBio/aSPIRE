# aSPIre
abundance of spliced peptides identified using rescoring

## overview and requirements
*aSPIre* provides the quantification of spliced and non-spliced peptides that were identified with *inSPIRE*.
In order to run *aSPIre*, you need access to a Windows machine that has [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view) installed. **Skyline** is a freely available tool for MS analysis.

*aSPIre* processes peptide-spectrum matches (PSMs) that were assigned by *inSPIRE*, quantifies them using Skyline and constructs a generation kinetic for each identified peptide.


## sample list
The user must provide *inSPIRE* final assignments and features. All information must be provided in the sample_list.csv (see an example below). You can edit the sample list using, e.g., MS Excel. In any case, make sure to save it as file with comma-separated values (.csv) and NOT as .xlsx notebook!
