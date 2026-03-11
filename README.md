# Tajima's D with K2P Modification

This project calculates **Tajima’s D** for DNA sequences in a FASTA format. It includes a standard calculation and a modified version using the **Kimura 2-Parameter (K2P)** model.

## Biological Realism Extension
To improve the biological realism of the standard Tajima's D statistic, I modified the calculation of average pairwise differences to use the K2P model instead of raw Hamming distance. 

Biologically, transition mutations occur more frequently than transversions. By applying the K2P model, we account for these differing mutation rates and estimate "invisible" substitutions that might be obscured by saturation over time.

## The script automates the following steps:

Data Acquisition: Queries the NCBI GenBank database (using Entrez) for specific organismal sequences (e.g., E. coli).

Population Simulation: Generates a synthetic population by applying stochastic SNP mutations to a base sequence, creating a controlled dataset for analysis.

Standard Tajima’s D Calculation: 

  * Parses the FASTA-formatted sequence alignment.Identifies segregating sites ($S$).
  * Calculates the mean number of pairwise differences using raw Hamming distance.
  * Computes the final $D$ statistic to detect signatures of selection or population expansion.

Biological Extension (K2P Model): Re-calculates the statistic using the Kimura 2-Parameter model to account for the biological reality that transitions and transversions occur at different rates.

## How to Run
1. Install dependencies: `pip install biopython numpy`
2. Run the script: `python tajima_analysis.py`
