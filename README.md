# Tajima's D with K2P Modification

This project calculates **Tajima’s D** for DNA sequences in a FASTA format. It includes a standard calculation and a modified version using the **Kimura 2-Parameter (K2P)** model.

## Biological Realism Extension
To improve the biological realism of the standard Tajima's D statistic, I modified the calculation of average pairwise differences to use the K2P model instead of raw Hamming distance. 

Biologically, transition mutations occur more frequently than transversions. By applying the K2P model, we account for these differing mutation rates and estimate "invisible" substitutions that might be obscured by saturation over time.

## How to Run
1. Install dependencies: `pip install biopython numpy`
2. Run the script: `python tajima_analysis.py`
